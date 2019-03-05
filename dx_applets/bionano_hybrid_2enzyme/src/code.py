import os
import glob
import dxpy

import dx_utils

BIONANO_ROOT = '/Solve3.2.1_04122018/'
SCRIPTS_DIR = os.path.join(BIONANO_ROOT, 'PIPELINE', 'Pipeline')

HYBRID_DIR = os.path.join(BIONANO_ROOT, 'HybridScaffold', '04122018')
TOOLS_DIR = os.path.join(BIONANO_ROOT, 'RefAligner', '7437.7523rel')


os.environ['PATH'] = HYBRID_DIR + \
    os.pathsep + os.environ['PATH']

# configure R library paths
os.link('/usr/lib/R/modules/lapack.so', '/usr/lib/R/lib/libRlapack.so')
os.link('/usr/lib/libblas.so', '/usr/lib/R/lib/libRblas.so')

@dxpy.entry_point("main")
def main(**job_inputs):
    bionano_cmap_1_link = job_inputs['bng_enzyme1']
    bionano_cmap_2_link = job_inputs['bng_enzyme2']
    ngs_fasta_link = job_inputs['ngs_fasta_or_cmap']
    args_xml_link = job_inputs.get('args_xml')

    # Download all the inputs
    bionano_cmap_1_filename = os.path.join(
        '/home/dnanexus', dx_utils.download_and_gunzip_file(bionano_cmap_1_link))
    bionano_cmap_2_filename = os.path.join(
        '/home/dnanexus', dx_utils.download_and_gunzip_file(bionano_cmap_2_link))
    ngs_fasta_filename = os.path.join('/home/dnanexus', dx_utils.download_and_gunzip_file(ngs_fasta_link))

    if args_xml_link:
        args_xml_filename = dx_utils.download_and_gunzip_file(args_xml_link)
    else:
        args_xml_filename = os.path.join(HYBRID_DIR, 'TGH', 'hybridScaffold_two_enzymes.xml')
    output_dir = "hybrid_scaffold_output"

    dx_utils.run_cmd(['mkdir', output_dir])
    results_tar = output_dir + '_results.tar'

    cmd = ["Rscript", os.path.join(HYBRID_DIR, "runTGH.R"), "--help"]
    dx_utils.run_cmd(cmd)

    scaffold_cmd = ["Rscript", os.path.join(HYBRID_DIR, "runTGH.R"), "-N", ngs_fasta_filename,
                    "-b1", bionano_cmap_1_filename, "-b2", bionano_cmap_2_filename, "-O", output_dir,
                    "-R", os.path.join(TOOLS_DIR, 'RefAligner'), "-t", results_tar,
                    "-e1", job_inputs['enzyme1_name'], "-e2", job_inputs['enzyme2_name']]

    if job_inputs.get("cuts1_file") and job_inputs.get("cuts2_file"):
        cuts1_file = dx_utils.download_and_gunzip_file(job_inputs["cuts1_file"])
        cuts2_file = dx_utils.download_and_gunzip_file(job_inputs["cuts2_file"])
        scaffold_cmd += ["-m1", cuts1_file, "-m2", cuts2_file]

    scaffold_cmd += [args_xml_filename]
    dx_utils.run_cmd(scaffold_cmd)

    scaffold_final = glob.glob(
        os.path.join(output_dir, 'TGH_M1', 'AGPExport', '*HYBRID*'))

    if not scaffold_final:
        print("ERROR: No hybrid scaffolds produced.")
        hybrid_scaffold_log = os.path.join(output_dir, 'TGH.log')
        dx_utils.run_cmd(["tail", "-n", "50", hybrid_scaffold_log])
    
    scaffold_final_ncbi = glob.glob(
        os.path.join(output_dir, 'hybrid_scaffolds*', '*_HYBRID_SCAFFOLD_NCBI.fasta'))[0]
    unscaffolded_final = glob.glob(
        os.path.join(output_dir, 'hybrid_scaffolds*', '*_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta'))[0]

    # make sure output files don't have colons
    dx_utils.run_cmd(["sed", "-i.bak", "s/:/_/g", scaffold_final_ncbi])
    dx_utils.run_cmd(["sed", "-i.bak", "s/:/_/g", unscaffolded_final])

    output = {
        "scaffold_fasta": [dxpy.dxlink(dxpy.upload_local_file(f)) for f in scaffold_final if f.endswith(".fasta")],
        "scaffold_output": [dxpy.dxlink(dxpy.upload_local_file(f)) for f in scaffold_final],
        "ncbi_scaffold_final": dx_utils.gzip_and_upload(scaffold_final_ncbi),
        "unscaffolded_final": dx_utils.gzip_and_upload(unscaffolded_final)
        }
    
    tar_name = "hybrid_scaffold_output.tar.gz"
    tar_cmd = "tar czvf {tar_name} {outdir}".format(
        tar_name=tar_name,
        outdir=output_dir)
    dx_utils.run_cmd(tar_cmd)
    output_id = dxpy.upload_local_file(tar_name)

    output["scaffold_targz"] = dxpy.dxlink(output_id)

    return output