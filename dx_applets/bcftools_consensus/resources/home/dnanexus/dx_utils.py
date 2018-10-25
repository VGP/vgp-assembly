from __future__ import print_function
import re
import os
import sys
import tempfile
import subprocess
import contextlib
import multiprocessing.pool

import dxpy


@contextlib.contextmanager
def set_env(**environ):
    """
    Temporarily set the process environment variables.

    >>> with set_env(PLUGINS_DIR=u'test/plugins'):
    ...   "PLUGINS_DIR" in os.environ
    True

    >>> "PLUGINS_DIR" in os.environ
    False

    :type environ: dict[str, unicode]
    :param environ: Environment variables to set
    """
    old_environ = dict(os.environ)
    os.environ.update(environ)
    try:
        yield
    finally:
        os.environ.clear()
        os.environ.update(old_environ)


def normalize_timedelta(timedelta):
    """
    Given a string like "1w" or "5d", convert it to an integer in seconds.
    Integers without a suffix are interpreted as seconds.
    Note: not related to the datetime timedelta class.
    """
    try:
        return int(timedelta)
    except ValueError as e:
        suffix_multipliers = {'s': 1, 'm': 60, 'h': 60*60, 'd': 60*60*24, 'w': 60*60*24*7,
                              'M': 60*60*24*30, 'y': 60*60*24*365}

        re_match = re.match('([0-9]+)([\s]*)([a-zA-Z]*)', timedelta)
        if re_match is None:
            msg = 'Could not parse time delta {0}'.format(timedelta)
            raise dxpy.AppInternalError(msg)
        
        t, space, suffix = re_match.groups()

        if suffix in suffix_multipliers:
            normalized_time = int(t) * suffix_multipliers[suffix]
        elif suffix not in suffix_multipliers and len(suffix) > 1 and suffix[0] in suffix_multipliers:
            old_suffix = suffix
            suffix = suffix[0]
            print('Not familiar with suffix {0}.  Assuming you meant {1}.'.format(old_suffix, suffix), file=sys.stderr)
            normalized_time = int(t) * suffix_multipliers[suffix]
        else:
            msg = 'Valid suffixes for time duration are {0}.'.format(str(suffix_multipliers.keys()))
            raise dxpy.AppInternalError(msg)

    return normalized_time


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
    if type(cmd) is list:
        shell=False
        executable = None
        print(subprocess.list2cmdline(cmd), file=sys.stderr)
    else:
        shell=True
        executable = '/bin/bash'
        print(cmd, file=sys.stderr)

    if returnOutput:
        output = subprocess.check_output(cmd, shell=shell, executable=executable).strip()
        print(output)
        return output
    else:
        subprocess.check_call(cmd, shell=shell, executable=executable)


def run_pipe(*cmds, **kwargs):
    """
    Function to run several commands that pipe to each other in a python 
    aware way.

    Valid cmd arguments are lists containing shell commands to execute.
    Valid kwargs are 'outputFile' and 'returnOutput'.

    Usage 1: Pipe multiple commands together and print output to file
        example_cmd1 = ['dx', 'download', 'file-xxxx']
        example_cmd2 = ['gunzip']
        out_f = "somefilename.fasta"
        run_pipe(example_cmd1, example_cmd2, outputFile=out_f)

        This function will print and execute the following command:
        'dx download file-xxxx | gunzip > somefilename.fasta'

    Usage 2: Pipe multiple commands together and return output
        example_cmd1 = ['gzip', 'file.txt']
        example_cmd2 = ['dx', 'upload', '-', '--brief']
        file_id = run_pipe(example_cmd1, example_cmd2, returnOutput=True)

        This function will print and execute the following command:
        'gzip file.txt | dx upload - --brief '
        and return the output.

    Usage 3: Pipe a single command with output to file
        example_cmd = ['ls', '.']
        out_f = 'ls.txt'
        run_pipe(example_cmd, outputFile=out_f)

        This function will print and execute the following command:
        'ls . > ls.txt'
    """
    # parse kwargs
    outputFile = kwargs.get('outputFile')
    returnOutput = kwargs.get('returnOutput', False)

    # copy provided commands into a list instead of a generator
    cmds = [c for c in cmds]

    # check that cmds are lists and not strings
    if type(cmds[0]) is not list:
        raise SyntaxError('Commands in pipe must be of type list')

    # pretty print the provided command
    cmd_str = _list2cmdlines_pipe(*cmds)
    if outputFile is not None:
        cmd_str += ' > {0}'.format(outputFile)
    print(cmd_str, file=sys.stderr)

    # there could be the case where there is only one command provided and we
    # want to stdout to file
    num_cmd = len(cmds)
    if num_cmd == 1 and outputFile is None:
        cmd = cmds[-1]
        if returnOutput is True:
            output = subprocess.check_output(cmd)
            return output.strip()
        else:
            subprocess.check_call(cmd)
            return
    elif num_cmd == 1 and outputFile:
        cmd = cmds[-1]
        with open(outputFile, "w") as fopen:
            subprocess.check_call(cmd, stdout=fopen)
        return

    # if multiple commands are provided, we need to pipe them together
    # intialize the first command
    cmd_process = []
    init_process = subprocess.Popen(cmds[0], stdout=subprocess.PIPE)
    cmd_process.append(init_process)
    # run all commands except the last one by piping
    for i in range(1, num_cmd - 1):
        # start the next command
        prev_process = cmd_process[i-1]
        process = subprocess.Popen(cmds[i],
                                   stdin=prev_process.stdout,
                                   stdout=subprocess.PIPE)
        # close the previous command -- this will take care of the case when
        # cmd_process[i] closes before cmd_process[-1] is done
        cmd_process.append(process)
        prev_process.stdout.close()

    # run the last command
    prev_process = cmd_process[-1]
    if outputFile is not None:
        with open(outputFile, "w") as fopen:
            subprocess.check_call(cmds[-1], stdin=prev_process.stdout, stdout=fopen)
            prev_process.stdout.close()
        return
    elif returnOutput is True:
        output = subprocess.check_output(cmds[-1], stdin=prev_process.stdout)
        prev_process.stdout.close()
        return output.strip()
    else:
        subprocess.check_call(cmds[-1], stdin=prev_process.stdout)
        prev_process.stdout.close()


def _list2cmdlines_pipe(*cmds):
    cmdline = ' | '.join([subprocess.list2cmdline(cmd) for cmd in cmds])

    return cmdline


class cd:
    '''From: http://stackoverflow.com/questions/431684/how-do-i-cd-in-python
       Context manager for changing the current working directory

       Usage:
       with cd():
           do_the_thing
           # this will create a temp directory with a randomly
           # generated name, doe the thing, then delete the temp dir

       with cd(my_file_dir):
           do_the_thing
           # this will do the thing in my_file_dir and not delete the directory

       with cd(tempDir=my_temp_dir):
           do_the_thing
           # this will create a temp dir with path my_temp_dir, do the thing,
           # then delete the temp dir
       '''
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
            subprocess.check_call(['rm', '-rf', self.newPath], shell=False)


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
    gzip_cmd = ['gzip', '--fast', '-c', fn]
    upload_cmd = ['dx', 'upload', '-', '--brief', '--path', rfn + '.gz']
    fid = run_pipe(gzip_cmd, upload_cmd, returnOutput=True)

    return dxpy.dxlink(fid)


def tar_files_and_upload(filenames, prefix):
    with tempfile.NamedTemporaryFile(delete=False) as fh:
        fh.write('\n'.join(filenames))

    ofn = '{0}.tar.gz'.format(prefix)
    tar_cmd = ['tar', 'cvf', '-', '--files-from', fh.name]
    gzip_cmd = ['gzip', '--fast']
    upload_cmd = ['dx', 'upload', '--brief', '--destination', ofn, '-']
    fid = run_pipe(tar_cmd, gzip_cmd, upload_cmd, returnOutput=True)
    cmd = ['rm', fh.name]
    run_cmd(cmd)

    return dxpy.dxlink(fid)


def remove_special_chars(string):
    '''function that replaces any characters in a string that are not
    alphanumeric or _ or . so that they do not cause issues in commands'''
    string = "".join(
        char for char in string if char.isalnum() or char in ['_', '.'])

    return string


def download_and_gunzip_file(input_file, skip_decompress=False, additional_pipe=None, create_named_pipe=False, input_filename=None):
    input_file = dxpy.DXFile(input_file)
    if input_filename is None:
        input_filename = input_file.describe()['name']
    ofn = remove_special_chars(input_filename)

    cmds = []
    cmds.append(['dx', 'download', input_file.get_id(), '-o', '-'])
    if input_filename.endswith('.tar.gz'):
        ofn = 'tar_output_{0}'.format(ofn.replace('.tar.gz', ''))
        cmds.append(['tar', '--no-same-owner', '-zxvf', '-'])
    elif input_filename.endswith('.tar.bz2'):
        ofn = 'tar_output_{0}'.format(ofn.replace('.tar.bz2', ''))
        cmds.append(['tar', '--no-same-owner', '-jxvf', '-'])
    elif input_filename.endswith('.tar'):
        ofn = 'tar_output_{0}'.format(ofn.replace('.tar', ''))
        cmds.append(['tar', '--no-same-owner', '-xvf', '-'])
    elif (os.path.splitext(input_filename)[-1] == '.gz') and not skip_decompress:
        cmds.append(['gunzip'])
        ofn = os.path.splitext(ofn)[0]
    if additional_pipe is not None:
        # check that additional pipe is not a string and does not have multiple
        # pipes in it
        if type(additional_pipe) is list:
            cmds.append(additional_pipe)
        else:
            # TODO: there could be an edge case where the additional pipe
            # contains the | character within quotation marks and it does not
            # signify pipe
            if additional_pipe.contains('|'):
                pipes = additional_pipe.split('|')
                pipe_cmds = [p.split() for p in pipes]
                cmds.extend(pipe_cmds)
            else:
                pipe = additional_pipe.split()
                cmds.append(pipe)

    if create_named_pipe:
        named_pipe_cmd = ['mkfifo', ofn]
        run_cmd(named_pipe_cmd)

    run_pipe(*cmds, outputFile=ofn)
    return ofn
