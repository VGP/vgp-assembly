import os
import subprocess
import glob
import dxpy


BIONANO_ROOT = '/Solve3.2.1_04122018/'
SCRIPTS_DIR = os.path.join(BIONANO_ROOT, 'PIPELINE', 'Pipeline')

HYBRID_DIR = os.path.join(BIONANO_ROOT, 'HybridScaffold', '04122018')
TOOLS_DIR = os.path.join(BIONANO_ROOT, 'RefAligner', '7437.7523rel')


os.environ['PATH'] = HYBRID_DIR + \
    os.pathsep + os.environ['PATH']


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
    bionano_cmap_link = job_inputs['refinefinal_merged_cmap']
    ngs_fasta_link = job_inputs['ngs_fasta_or_cmap']
    args_xml_link = job_inputs.get('args_xml')

    # Download all the inputs
    bionano_cmap_filename = download_and_gunzip_file(bionano_cmap_link)
    ngs_fasta_filename = download_and_gunzip_file(ngs_fasta_link)
    if args_xml_link:
        args_xml_filename = download_and_gunzip_file(args_xml_link)
    else:
        args_xml_filename = os.path.join(HYBRID_DIR, 'hybridScaffold_config.xml')

    output_dir = "hybrid_scaffold_output"

    scaffold_cmd = (
        "perl {dir}/hybridScaffold.pl -n {ngs_fasta} -b {cmap} "
        "-o {outdir} -c {args_xml} "
        .format(dir=HYBRID_DIR, ngs_fasta=ngs_fasta_filename, cmap=bionano_cmap_filename,
                args_xml=args_xml_filename, outdir=output_dir))
    scaffold_cmd += '-r {refaligner} '.format(refaligner=os.path.join(TOOLS_DIR, 'RefAligner'))

    if "conflict_resolution_file" in job_inputs:
        conflict_resolution_file = download_and_gunzip_file(job_inputs["conflict_resolution_file"])
        scaffold_cmd += '-M {cr_file} '.format(conflict_resolution_file)
    else:
        scaffold_cmd += '-B {b_level} -N {n_level} '.format(
            b_level=job_inputs["b_conflict_filter"], n_level=job_inputs["n_conflict_filter"])

    if job_inputs["generate_molecules"] is True:
        scaffold_cmd += '-x '
        scaffold_cmd += '-p {0}'.format(SCRIPTS_DIR)

        try:
            molecules_bnx_file = download_and_gunzip_file(job_inputs["molecules_bnx_file"])
            scaffold_cmd += '-m {0} '.format(molecules_bnx_file)

        except KeyError:
            raise dxpy.AppError("Molecules BNX file required for Align Molecules flag (-x)")

        try:
            optargs_xml = download_and_gunzip_file(job_inputs["optargs_xml"])
            scaffold_cmd += '-q {0} '.format(optargs_xml)

        except KeyError:
            raise dxpy.AppError("OptArgs XML file required for Align Molecules flag (-x)")

    if job_inputs["generate_chimeric"] is True:
        scaffold_cmd += '-y '

        if molecules_bnx_file:
            scaffold_cmd += '-m {0} '.format(molecules_bnx_file)

        else:
            try:
                molecules_bnx_file = download_and_gunzip_file(job_inputs["molecules_bnx_file"])
                scaffold_cmd += '-m {0} '.format(molecules_bnx_file)

            except KeyError:
                raise dxpy.AppError("Molecules BNX file required for Generate Molecules flag")

        if "err_files" in job_inputs:
            err_files = [download_and_gunzip_file(err_file) for err_file in job_inputs["err_files"]]
            err_cmd = ' '.join(['-e {0}'.format(err) for err in err_files])
            scaffold_cmd += err_cmd
    run_cmd(scaffold_cmd)

    run_cmd('tree {0}'.format(output_dir))

    tar_name = "hybrid_scaffold_output.tar.gz"
    tar_cmd = "tar czvf {tar_name} {outdir}".format(
        tar_name=tar_name,
        outdir=output_dir)
    run_cmd(tar_cmd)
    output_id = dxpy.upload_local_file(tar_name)

    scaffold_final = glob.glob(
        os.path.join(output_dir, 'hybrid_scaffolds*', '*_HYBRID_SCAFFOLD.fasta'))
    scaffold_final.extend(glob.glob(
        os.path.join(output_dir, 'hybrid_scaffolds*', '*_HYBRID_SCAFFOLD.cmap')))
    scaffold_final.extend(glob.glob(
        os.path.join(output_dir, 'hybrid_scaffolds*', '*_HYBRID_SCAFFOLD.agp')))
    scaffold_output = glob.glob(
        os.path.join(output_dir, 'hybrid_scaffolds*', '*_HYBRID_SCAFFOLD*'))
    cut_and_conflict = glob.glob(os.path.join(output_dir, 'hybrid_scaffolds*', 'conflicts*.txt'))
    cut_and_conflict.extend(glob.glob(os.path.join(output_dir, 'hybrid_scaffolds*', '*_annotations.bed')))
    return {"scaffold_final": [dxpy.dxlink(dxpy.upload_local_file(f)) for f in scaffold_final],
            "scaffold_output": [dxpy.dxlink(dxpy.upload_local_file(f)) for f in scaffold_output if f not in scaffold_final],
            "scaffold_targz": dxpy.dxlink(output_id),
            "cut_and_conflict": [dxpy.dxlink(dxpy.upload_local_file(f)) for f in cut_and_conflict]}
