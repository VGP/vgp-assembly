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
    Context manager generator to temporarily set the subprocess environment variables.

    Args:
        environ (dict): Environment variable to set

    Yields:
        An environment with environment variables set as specified.
        On exit, the environment will return to previous configuration.

    Examples:
        Usage 1: Set environment variable
        # inside environment
        >>> with set_env(PLUGINS_DIR=u'test/plugins'):
        ...    "PLUGINS_DIR" in os.environ
        True

        # outside environment
        >>> "PLUGINS_DIR" in os.environ
        False

        Usage 2: Unset environment variable
        >>> with set_env(PYTHONPATH=''):
        ...    print(os.environ["PYTHONPATH"])
        <BLANKLINE>

        Usage 3: Manipulate multiple variables
        >>> myenv = {"PLUGINS_DIR": u'test/plugins', "PYTHONPATH": u'some/python/path'}
        >>> with set_env(**myenv):
        ...   print(os.environ["PLUGINS_DIR"])
        ...   print(os.environ["PYTHONPATH"])
        test/plugins
        some/python/path
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
    Integers without a suffix are interpreted as seconds. It is assumed this
    function is run inside a dnanexus job.

    Note: not related to the datetime timedelta class.

    Args:
        timedelta (int or str): An integer in seconds or a string of the form "1w" or "5d"

    Returns:
        normalized_time (int): An integer in seconds

    Raises:
        dxpy.AppInternalError if suffixes are not known. Known suffixes are
        defined in the dictionary suffix_multipliers.

    Examples:
        >>> normalize_timedelta(500)
        500

        >>> normalize_timedelta('5w')
        3024000

        >>> normalize_timedelta('1d')
        86400

        >>> normalize_timedelta('1h')
        3600

    """
    try:
        return int(timedelta)
    except ValueError as e:
        suffix_multipliers = {'s': 1, 'm': 60, 'h': 60*60, 'd': 60*60*24, 
                              'w': 60*60*24*7,
                              'M': 60*60*24*30, 'y': 60*60*24*365}

        re_match = re.match('([0-9]+)([\s]*)([a-zA-Z]*)', timedelta)
        if re_match is None:
            msg = 'Could not parse time delta {0}'.format(timedelta)
            raise dxpy.AppInternalError(msg)

        t, _, suffix = re_match.groups()

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
    """Internal function used in the class ExceptionAwarePool."""
    try:
        result = func(*args, **kwargs)
    except Exception as e:
        q.put(e.__str__())

    return result

def get_dxlink_filesizes(dx_links):
    """Run dx describe on a list of DNAnexus dxlink inputs to get the
    corresponding file sizes.
    
    Args:
        dx_links (list of dicts): dxlink dicts containing '$dnanexus_link' as 
            key and file-id as value
    
    Returns:
        list: corresponding filesizes in bytes, output of 'dx describe'
        command
    """
    input = {'objects': [file['$dnanexus_link'] for file in dx_links]}
    descriptions = dxpy.api.system_describe_data_objects(input, always_retry=True)

    sizes = [d['describe']['size'] for d in descriptions['results']]

    return sizes

def get_project(project_name):
    """Try to find the project with the given name or id on DNAnexus.
    It is assumed that the user or job is logged in prior to running this 
    function so dx API queries are returned."""

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
    """
    A wrapper function around multiprocessing Pool class to catch exceptions
    that occur internal to the pool.

    Args:
        Same *args and **kwargs as multiprocessing.pool.Pool

    Functions:
        EAP.apply_async() and EAP.join() act the same way as Pool.apply_async()
        and Pool.apply_join() but will raise an exception if one is encountered
        in any of the pooled functions.
    """
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
    """
    Queries a worker's /proc/meminfo for available memory and returns a float
    of the specified suffix size

    Args:
        suffix (str): One of 'M', 'K' or 'G' to return memory in Mib, KiB or
                      GiB, respectively.

    Returns:
        float: total_memory read from meminfo in MiB, KiB or GiB
            depending on specified suffix.
    Raises:
        dxpy.DXError is raised if suffix is not recognized or system memory
        cannot be read.
    """
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


def run_cmd(cmd, return_output=False, output_file=None):
    """Print the cmd to stderr and execute using a subprocess.

    Args:
        cmd (str or list): Command to be executed using subprocess. Input of 
            type 'list' is recommended. When input is of type 'string',
            command is executed using /bin/bash and shell=True is specified.
        return_output (bool): If command output should be returned.
        output_file (str): Filename of the file command std output should be redirected to.

    Note:
        return_output and output_file inputs are mutually exclusive.
        If output_file exists, it will be overwritten.

    Returns:
        Before executing, 'cmd' input is pretty-printed to stderr.
        If return_output is provided, the output string is returned.
        If output_file is provided, stdout is redirected to the specified file
        and the function returns None.

    Raises:
        CalledProcessError may be raised if the command failed to execute.

    Examples:
        run_cmd(['echo', 'hello world'], output_file='test1.txt')

    """
    if type(cmd) is list:
        shell = False
        executable = None
        if output_file:
            print_cmd = cmd + ['>', output_file]
        else:
            print_cmd = cmd
        print(subprocess.list2cmdline(print_cmd), file=sys.stderr)
    else:
        shell = True
        executable = '/bin/bash'
        if output_file:
            print_cmd = cmd + '>' + output_file
        else:
            print_cmd = cmd
        print(print_cmd, file=sys.stderr)

    if return_output:
        output = subprocess.check_output(cmd, shell=shell, executable=executable).strip()
        print(output)
        return output.strip()
    elif output_file:
        with open(output_file, 'w') as fopen:
            subprocess.Popen(cmd, shell=shell, executable=executable, stdout=fopen)
    else:
        subprocess.check_call(cmd, shell=shell, executable=executable)


def run_pipe(*cmds, **kwargs):
    """
    Function to run several commands that pipe to each other in a python 
    aware way.

    Args:
        *cmds: any number of inputs of cmds (type list) to pipe together
        **kwargs: The only valid kwargs are output_file and return_output. Other
        arbitrary kwargs are ignored.
            output_file (str): filename to redirected stdout to at end of pipe
            return_output (bool): Whether to return stdout at end of pipe

    Note: 
        output_file and return_output kwargs are mutually exclusive.
        If output_file exists, it will be overwritten.

    Returns:
        Piped command is pretty-printed to stderr.
        output (str) is returned if return_output=True is passed as **kwarg.
        If output_file is passed as kwarg, stdout is redirected to the given file.

    Raises:
        TypeError: if *cmds specified are not of type list
        subprocess.CalledProcessError: if any subprocess in pipe returns exit
            code not 0.

    Examples:
        Usage 1: Pipe multiple commands together and print output to file
            example_cmd1 = ['dx', 'download', 'file-xxxx']
            example_cmd2 = ['gunzip']
            out_f = "somefilename.fasta"
            run_pipe(example_cmd1, example_cmd2, output_file=out_f)

            This function will print and execute the following command:
            'dx download file-xxxx | gunzip > somefilename.fasta'

        Usage 2: Pipe multiple commands together and return output
            example_cmd1 = ['gzip', 'file.txt']
            example_cmd2 = ['dx', 'upload', '-', '--brief']
            file_id = run_pipe(example_cmd1, example_cmd2, return_output=True)

            This function will print and execute the following command:
            'gzip file.txt | dx upload - --brief '
            and return the output.

        Usage 3: Pipe a single command with output to file
            run_pipe(['echo', 'hello world'], output_file='test2.txt')
            Note: This calls the run_cmd function instead of run_pipe.

        Usage 4: A command failing mid-pipe should return CalledProcessedError
            >>> run_pipe(['echo', 'hi:bye'], ['grep', 'blah'], ['cut', '-d', ':', '-f', '1'])
            Traceback (most recent call last):
                  ...
            CalledProcessError: Command '['grep', 'blah']' returned non-zero exit status 1
    """
    # parse kwargs
    output_file = kwargs.get('output_file')
    return_output = kwargs.get('return_output', False)

    # copy provided commands into a list instead of a generator
    cmds = [c for c in cmds]

    # check that cmds are lists and not strings
    if type(cmds[0]) is not list:
        raise TypeError('Commands in pipe must be of type list')

    num_cmd = len(cmds)
    # if only one command is provided, use run_cmd instead
    if num_cmd == 1:
        return run_cmd(cmds[0], return_output, output_file)

    # from here we can assume multiple commands are provided
    # pretty print the provided command
    cmd_str = list2cmdlines_pipe(*cmds)
    if output_file is not None:
        cmd_str += ' > {0}'.format(output_file)
    print(cmd_str, file=sys.stderr)

    # if multiple commands are provided, we need to pipe them together
    # initialize the first command
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
    output = None
    prev_process = cmd_process[-1]
    if output_file is not None:
        with open(output_file, "w") as fopen:
            subprocess.check_call(cmds[-1], stdin=prev_process.stdout, stdout=fopen)
            prev_process.stdout.close()
    elif return_output is True:
        output = subprocess.check_output(cmds[-1], stdin=prev_process.stdout)
        prev_process.stdout.close()
        output = output.strip()
    else:
        subprocess.check_call(cmds[-1], stdin=prev_process.stdout)
        prev_process.stdout.close()

    # check that all intermediate commands finished successfully
    for i in range(len(cmd_process)):
        cmd = cmds[i]
        curr_proc = cmd_process[i]
        # Polling is needed first in order to set the returncode attr
        curr_proc.poll()
        returncode = curr_proc.returncode
        if returncode != 0 and returncode is not None:
            raise subprocess.CalledProcessError(
                returncode,  cmd, output=None)

    # if return_output is True, then the variable output contains the subprocess
    # output. If return_output is False, the variable output should be set to 'None',
    # which is the same behavior as a simple `return` statement.
    return output


def list2cmdlines_pipe(*cmds):
    """Given a list of cmds (list of lists), pretty-print the commands as a shell
    pipe string.

    Args:
        *cmds: any number of inputs of cmds (type list) to pipe together

    Returns:
        string: command lists converted to string and joined using ' | '.

    Examples:
        >>> list2cmdlines_pipe(['cat', 'file.txt'], ['gzip'], ['dx', 'upload', '-'])
        'cat file.txt | gzip | dx upload -'
    """
    cmdline = ' | '.join([subprocess.list2cmdline(cmd) for cmd in cmds])

    return cmdline


class cd:
    '''
    Context manager for changing the current working directory

    Args:
        new_path (string): Optional, specify path to cd to
        temp_dir (string): Optional, specify temporary directory to create and
            cd to

    Note:
        If no args specified, cd() will create an arbitary temp dir and cd to it
        If arg is specified without a keyword, it will be assumed as 'new_path'.

    Yields:
        Upon entry, context will be set to the specified directory.
        Upon exit, directory specified in temp_dir or directory created when no
        args are specified is deleted. If new_path is specified, it is not deleted.

    Source: http://stackoverflow.com/questions/431684/how-do-i-cd-in-python

    Examples:
       with cd():
           do_the_thing
           # this will create a temp directory with a randomly
           # generated name, doe the thing, then delete the temp dir

       with cd(my_file_dir):
           do_the_thing
           # this will do the thing in my_file_dir and not delete the directory

       with cd(temp_dir=my_temp_dir):
           do_the_thing
           # this will create a temp dir with path my_temp_dir, do the thing,
           # then delete the temp dir
       '''
    def __init__(self, new_path=None, temp_dir=None):
        if new_path is not None:
            self.newPath = new_path
            self.removeFolder = False
        else:
            self.newPath = tempfile.mkdtemp(dir=temp_dir)
            self.removeFolder = True

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
        if self.removeFolder:
            subprocess.check_call(['rm', '-rf', self.newPath], shell=False)


def schedule_lpt(jobs, num_bins):
    """This function implements the Longest Processing Time algorithm to get
    a good division of labor for the multiprocessor scheduling problem.

    Args:
        jobs (dict or list): A dictionary with string 'key' specifying job
            name and float 'value' specifying job weight (how long job should 
            run compared to other jobs). A list of tuples may also be passed 
            if it is equivalent to dict.items()
        num_bins (int): Number of groups to split jobs into.

    Returns:
        List of lists: Each group (list) in list is a group of jobs that should
            be run together on a single worker or instance and consists of 
            tuples as provided in jobs input.

    Examples:
        # it's assumed that there's an app-specific way to generate a filelist
        # of filenames and sizes
        fl = filenames_and_sizes(files)
        fl_groups = schedule_lpt(fl, num_jobs)
        for group in fl_groups:
            print group
            job = dxpy.new_dxjob({'files': group}, 'subjob_name')
            output['output_files'].append(job.get_output_ref('output_files'))
    """

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


def gzip_and_upload(fn, rfn=None, compression_level=1):
    """This function is shorthand for running 'gzip' and 'dx upload' using
    subprocess on a given local file.

    Args:
        fn (string) = local filename or location
        rfn (string) = Optional, output filename to upload fn to.
        compression_level (int) = Level of compression between 1 and 9 to compress
            file to. Specify 1 for gzip --fast and 9 for gzip --best. If not 
            specified, --fast is assumed.

    Returns:
        dict: DNAnexus link pointing to the uploaded file

    Raises:
        ValueError: if compression_level not between 1 and 9
        CalledProcessError: propogated from run_pipe if called command fails
    
    Examples:
        >>> # create an example file
        >>> open("testfile1.txt", "w").write("my test file")
        >>> # gzip and upload
        >>> fid = gzip_and_upload("testfile1.txt", compression_level=6)
        >>> # confirm that fid is a dxlink dict
        >>> type(fid) is dict
        True
        >>> "$dnanexus_link" in fid
        True
    """
    if rfn is None:
        rfn = os.path.split(fn)[-1]

    # check if output filename ends with .gz, if not, add it
    if not rfn.endswith('.gz'):
        rfn += '.gz'

    # check that compression level is valid
    if compression_level < 1 or compression_level > 9:
        raise ValueError('Compression level must be between 1 and 9')

    gzip_cmd = ['gzip', '-{0}'.format(compression_level), '-c', fn]
    upload_cmd = ['dx', 'upload', '-', '--brief', '--path', rfn]
    fid = run_pipe(gzip_cmd, upload_cmd, returnOutput=True)

    return dxpy.dxlink(fid)


def tar_files_and_upload(filenames, prefix, compression_level=1):
    """This function is shorthand for running 'tar', 'gzip', and 'dx upload'
    using subprocess on a list of given local files.

    Args:
        filenames (list of strings) = local filenames or locations
        prefix (string) = name to give to output tar archive
        compression_level (int) = Level of compression between 1 and 9 to 
            compress tar to. Specify 1 for gzip --fast and 9 for gzip --best.
            If not specified, --fast is assumed.

    Returns:
        dict: DNAnexus link pointing to the uploaded tar archive

    Raises:
        ValueError: if compression_level not between 1 and 9
        CalledProcessError: propogated from run_pipe if called command fails
    
    Examples:
        >>> # create a couple test files
        >>> open("test_file1.txt", "w").write("testfile1")
        >>> open("test_file2.txt", "w").write("testfile2")
        >>> # tar the files and upload
        >>> fid = tar_files_and_upload(["test_file1.txt", "test_file2.txt"], "mytesttar")
        >>> # confirm that fid is a dict with "$dnanexus_link" key
        >>> type(fid) is dict
        True
        >>> "$dnanexus_link" in fid
        True
    """
    # create a temporary file as list of filenames
    with tempfile.NamedTemporaryFile(delete=False) as fh:
        fh.write('\n'.join(filenames))

    # check that compression level is valid
    if compression_level < 1 or compression_level > 9:
        raise ValueError('Compression level must be between 1 and 9')

    # specify output tar archive with provided prefix
    ofn = '{0}.tar.gz'.format(prefix)

    # generate commands to pipe
    tar_cmd = ['tar', 'cvf', '-', '--files-from', fh.name]
    gzip_cmd = ['gzip', '-{0}'.format(compression_level)]
    upload_cmd = ['dx', 'upload', '--brief', '--destination', ofn, '-']

    # run pipe
    fid = run_pipe(tar_cmd, gzip_cmd, upload_cmd, returnOutput=True)

    # remove temporary filelist file
    cmd = ['rm', fh.name]
    run_cmd(cmd)

    # return the DNAnexus link
    return dxpy.dxlink(fid)


def remove_special_chars(string):
    '''function that replaces any characters in a string that are not
    alphanumeric or _ or . so that they do not cause issues in commands

    Args:
        string (str): Any string with or without special characters

    Returns:
        string: Same as input but with special characters removed.'''
    string = "".join(
        char for char in string if char.isalnum() or char in ['_', '.'])

    return string


def download_and_gunzip_file(input_file, input_filename=None,
                             skip_decompress=False, additional_pipe=None,
                             create_named_pipe=False):
    """
    Shorthand for running dx download on a given input_file dx file link.
    Will additionally use subprocess to decompress or untar the file 
    automatically based on the name suffix of the file provided.

    Args:
        input_file (dict or string): DNAnexus link or file-id of file to download
        input_filename (string): Local filename to save file to. If not provided,
            platform filename is used.
        skip_decompress (bool): Whether to skip decompressing files of type *.gz
        additional_pipe (list or string): Additional commands to pipe the download
            command to. List input is recommended, string should include
            commands separated by the character '|'.
        create_named_pipe (bool): Whether to convert input_filename into a named
            pipe instead of a local output file.

    Notes:
        Supported filetypes: *.tar.gz, *.tar, *.tar.bz2, *.gz
        The arg skip_decompress only impacts files of type *.gz
        List input for additional_pipe is highly recommended over string.

    Returns:
        string: filename or local named pipe which file was downloaded to

    #TODO:
        - This function is kind of overloaded
        - Having '|' in additional_pipe which does not mean pipe can create a
        weird edge case which is not accounted for
    """
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

    run_pipe(*cmds, output_file=ofn)
    return ofn


if __name__ == "__main__":
    import doctest
    test_failures = doctest.testmod()[0]
    if test_failures > 0:
        print("Encountered {0} failures".format(test_failures))
        sys.exit(1)
    else:
        print("All tests passed.")