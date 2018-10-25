import os
import subprocess
import glob
import dxpy
import dx_utils

BIONANO_ROOT = '/Solve3.2.1_04122018/'
SCRIPTS_DIR = os.path.join(BIONANO_ROOT, 'PIPELINE', 'Pipeline')

HYBRID_DIR = os.path.join(BIONANO_ROOT, 'HybridScaffold', '04122018')
TOOLS_DIR = os.path.join(BIONANO_ROOT, 'RefAligner', '7437.7523rel')

os.environ['PATH'] = HYBRID_DIR + \
    os.pathsep + os.environ['PATH']


# ------------------------------------------------------
# /Solve_03062017Rel/HybridScaffold/03062017/scripts/fa2cmap_multi_color.pl

# Usage: /usr/bin/perl fa2cmap_multi_color.pl [options] <Args>
# Options:
#   -h : This help message
#   -v : Print program version information
#   -i : Input fasta file (Required)
#   -k : Input key file (If provided, use as contigIDs rather than sequentially)
#   -o : Output folder (Default: the same as the input file)
#   -e : Name or sequence of the enzyme (or a combination of name and sequence
#        in the format of name:sequence) followed by channel # (Can be multiple)
#   -m : Filter: Minimum nicks (Integer, default: 0)
#   -M : Filter: Minimum size (kb) (Integer, default: 0)
#   -B : Output BED file of Nbase gaps (Default: OFF)
#   -g : Minimum N gap size (bp) when generating Nbase gap file (Default: 1000)
#   -W : For web use only, and must be the first option (Default: OFF)
#   -S : Add an additional column, Strand, to the cmap file (Default: OFF)
#
# NOTE: CMAP index is 1-based, and is color-aware.
# -------------------------------------------------------

def run_cmd(cmd, returnOutput=False):
    print cmd
    if returnOutput:
        output = subprocess.check_output(
            cmd, shell=True, executable='/bin/bash').strip()
        print output
        return output
    else:
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')


def remove_special_chars(string):
    '''function that replaces any characters in a string that are not alphanumeric or _ or .'''
    string = "".join(
        char for char in string if char.isalnum() or char in ['_', '.'])

    return string


def download_and_gunzip_file(input_file, skip_decompress=False, additional_pipe=None):
    input_file = dxpy.DXFile(input_file)
    input_filename = input_file.describe()['name']
    ofn = remove_special_chars(input_filename)

    cmd = 'dx download ' + input_file.get_id() + ' -o - '
    if input_filename.endswith('.tar.gz'):
        ofn = 'tar_output_{0}'.format(ofn.replace('.tar.gz', ''))
        cmd += '| tar -zxvf - '
    elif (os.path.splitext(input_filename)[-1] == '.gz') and not skip_decompress:
        cmd += '| gunzip '
        ofn = os.path.splitext(ofn)[0]
    if additional_pipe is not None:
        cmd += '| ' + additional_pipe
    cmd += ' > ' + ofn
    print cmd
    subprocess.check_call(cmd, shell=True)

    return ofn


@dxpy.entry_point("main")
def main(**job_inputs):
    ngs_fasta_link = job_inputs['ngs_fasta']

    # Download all the inputs
    ngs_fasta_filename = download_and_gunzip_file(ngs_fasta_link)

    help_cmd = ('perl {0}/scripts/fa2cmap_multi_color.pl -v'.format(HYBRID_DIR))
    run_cmd(help_cmd)

    fa2cmap_cmd = (
        "perl {dir}/scripts/fa2cmap_multi_color.pl -i {ngs_fasta} -e {enzyme} {channel}"
        .format(dir=HYBRID_DIR, ngs_fasta=ngs_fasta_filename, 
                enzyme=job_inputs['enzyme_name'], 
                channel=job_inputs['channel_num'])
                   )

    prefix = ngs_fasta_filename.replace('.fasta', '').replace('.fa', '')
    print(prefix)
    dx_utils.run_cmd(fa2cmap_cmd)
    dx_utils.run_cmd('ls .')

    output = glob.glob(prefix + '*.cmap')[0]
    return {"reference_cmap": dxpy.dxlink(dxpy.upload_local_file(output))}
