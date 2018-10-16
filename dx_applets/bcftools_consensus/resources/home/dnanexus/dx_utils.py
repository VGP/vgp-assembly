from __future__ import print_function
import re
import os
import sys
import tempfile
import subprocess
import multiprocessing
import multiprocessing.pool

import dxpy

def normalize_timedelta(timedelta):
    """
    Given a string like "1w" or "5d", convert it to an integer in seconds.
    Integers without a suffix are interpreted as seconds.
    Note: not related to the datetime timedelta class.
    """
    try:
        return int(timedelta)
    except ValueError as e:
        t, suffix = timedelta[:-1], timedelta[-1:]
        suffix_multipliers = {'s': 1, 'm': 60, 'h': 60*60, 'd': 60*60*24, 'w': 60*60*24*7,
                              'M': 60*60*24*30, 'y': 60*60*24*365}
        if suffix not in suffix_multipliers:
            msg = 'Valid suffixes for time duration are shdwMy'
            raise dxpy.AppInternalError(msg)

        return int(t) * suffix_multipliers[suffix]


def _eap_wrapper(func, q, *args, **kwargs):
    try:
        result = func(*args, **kwargs)
    except Exception as e:
        q.put(e.__str__())
        pass

    return result


def get_project(project_name):
    '''Try to find the project with the given name or id.'''

    # First, see if the project is a project-id.
    try:
        project = dxpy.DXProject(project_name)
        return project
    except dxpy.DXError:
        pass

    project = dxpy.find_projects(name=project_name, name_mode='glob', return_handler=True, level="VIEW")
    project = [p for p in project]
    if len(project) < 1:
        print('Did not find project {0}.'.format(project_name), file=sys.stderr)
        sys.exit(1)
    elif len(project) > 1:
        print('Found more than 1 project matching {0}.'.format(project_name), file=sys.stderr)
        sys.exit(1)
    else:
        project = project[0]

    return project


class ExceptionAwarePool(multiprocessing.pool.Pool):
    def __init__(self, *args, **kwargs):
        super(ExceptionAwarePool, self).__init__(*args, **kwargs)
        self.manager = multiprocessing.Manager()
        self.q = self.manager.Queue()

    def apply_async(self, func, args=(), kwds={}, callback=None):
        return super(ExceptionAwarePool, self).apply_async(_eap_wrapper, args=(func, self.q) + args, kwds=kwds, callback=callback)

    def join(self):
        super(ExceptionAwarePool, self).join()
        if self.q.qsize() > 0:
            raise Exception('\n'.join([self.q.get() for i in xrange(self.q.qsize())]))


def get_memory(suffix='M'):
    if suffix == 'K':
        shift = 1
    elif suffix == 'M':
        shift = 1 << 10
    elif suffix == 'G':
        shift = 1 << 20
    else:
        raise dxpy.DXError('Unknown memory suffix {0}.  Please choose from K, M, or G.'.format(suffix))

    # Calc amount of memory available for gatk and Picard.
    total_mem = re.findall('^MemTotal:[\s]*([0-9]*) kB',
                           open('/proc/meminfo').read())
    if(len(total_mem) != 1):
        raise dxpy.DXError('Problem reading system memory from /proc/meminfo')
    return float(total_mem[0]) / shift


def run_cmd(cmd, returnOutput=False):
    print(cmd)
    if returnOutput:
        output = subprocess.check_output(cmd, shell=True, executable='/bin/bash').strip()
        print(output)
        return output
    else:
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')


class cd:
    '''From: http://stackoverflow.com/questions/431684/how-do-i-cd-in-python
       Context manager for changing the current working directory'''
    def __init__(self, newPath=None, tempDir=None):
        if newPath is not None:
            self.newPath = newPath
            self.removeFolder = False
        else:
            self.newPath = tempfile.mkdtemp(dir=tempDir)
            self.removeFolder = True

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
        if self.removeFolder:
            subprocess.check_call('rm -rf {0}'.format(self.newPath), shell=True)


def schedule_lpt(jobs, num_bins):
    '''This function implements the Longest Processing Time algorithm to get
    a good division of labor for the multiprocessor scheduling problem.'''

    def _index_min(values):
        # Efficient index of min from stackoverflow.com/questions/2474015
        return min(xrange(len(values)), key=values.__getitem__)

    # We expect a list of tuples, with the first value the name of the
    # job and the second value the weight.  If we are given a dict
    # then convert keys to job names and values to weights.
    if(type(jobs) == dict):
        jobs = zip(jobs.keys(), jobs.values())

    num_bins = min(num_bins, len(jobs))
    jobs.sort(key=lambda j: j[1], reverse=True)
    partition = {'groups': [[] for i in xrange(num_bins)],
                 'size': [0 for i in xrange(num_bins)]}

    for job in jobs:
        idx = _index_min(partition['size'])
        partition['groups'][idx] += [job[0]]
        partition['size'][idx] += job[1]

    return partition['groups']


def gzip_and_upload(fn, rfn=None):
    if rfn is None:
        rfn = os.path.split(fn)[-1]
    cmd = 'gzip --fast -c {0} | dx upload --brief --path {1}.gz -'.format(fn, rfn)
    print(cmd)
    fid = subprocess.check_output(cmd, shell=True).strip()

    return dxpy.dxlink(fid)


def tar_files_and_upload(filenames, prefix):
    with tempfile.NamedTemporaryFile(delete=False) as fh:
        fh.write('\n'.join(filenames))

    ofn = '{0}.tar.gz'.format(prefix)
    cmd = 'tar cvf - --files-from {0} | gzip --fast | dx upload --brief --destination {1} - '.format(fh.name, ofn)
    fid = run_cmd(cmd, returnOutput=True)
    cmd = 'rm {0} '.format(fh.name)
    run_cmd(cmd)

    return dxpy.dxlink(fid)


def download_and_gunzip_file(input_file, skip_decompress=False, additional_pipe=None, create_named_pipe=False, input_filename=None):
    input_file = dxpy.DXFile(input_file)
    if input_filename is None:
        input_filename = input_file.describe()['name']
    ofn = input_filename

    cmd = 'dx download ' + input_file.get_id() + ' -o - '
    if input_filename.endswith('.tar.gz'):
        ofn = 'tar_output_{0}'.format(ofn.replace('.tar.gz', ''))
        cmd += '| tar -zxvf - '
    elif (os.path.splitext(input_filename)[-1] == '.gz') and not skip_decompress:
        cmd += '| gunzip '
        ofn = os.path.splitext(ofn)[0]
    if additional_pipe is not None:
        cmd += '| ' + additional_pipe
    cmd += ' > "{0}"'.format(ofn)

    if create_named_pipe:
        named_pipe_cmd = 'mkfifo {0}'.format(ofn)
        run_cmd(named_pipe_cmd)
        cmd += '&'

    run_cmd(cmd)

    return ofn
