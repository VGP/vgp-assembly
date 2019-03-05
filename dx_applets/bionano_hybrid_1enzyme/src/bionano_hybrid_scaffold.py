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


@dxpy.entry_point("main")
def main(**job_inputs):
    bionano_cmap_link = job_inputs['refinefinal_merged_cmap']
    ngs_fasta_link = job_inputs['ngs_fasta_or_cmap']
    args_xml_link = job_inputs.get('args_xml')

    # Download all the inputs
    bionano_cmap_filename = dx_utils.download_and_gunzip_file(bionano_cmap_link)
    ngs_fasta_filename = dx_utils.download_and_gunzip_file(ngs_fasta_link)
    if args_xml_link:
        args_xml_filename = dx_utils.download_and_gunzip_file(args_xml_link)
    else:
        args_xml_filename = os.path.join(HYBRID_DIR, 'hybridScaffold_config.xml')

    output_dir = "hybrid_scaffold_output"

    scaffold_cmd = ["perl", os.path.join(HYBRID_DIR, "hybridScaffold.pl"), "-n", ngs_fasta_filename,
                    "-b", bionano_cmap_filename, "-o", output_dir, "-c", args_xml_filename,
                    "-r", os.path.join(TOOLS_DIR, 'RefAligner')]

    if "conflict_resolution_file" in job_inputs:
        conflict_resolution_file = dx_utils.download_and_gunzip_file(job_inputs["conflict_resolution_file"])
        scaffold_cmd += ["-M", conflict_resolution_file]
    else:
        scaffold_cmd += ["-B", job_inputs["b_conflict_filter"], "-N", job_inputs["n_conflict_filter"]]

    molecules_bnx_file = None
    if job_inputs["generate_molecules"] is True:
        scaffold_cmd += ["-x", "-p", SCRIPTS_DIR]

        try:
            molecules_bnx_file = dx_utils.download_and_gunzip_file(job_inputs["molecules_bnx_file"])
            scaffold_cmd += ["-m", molecules_bnx_file]

        except KeyError:
            raise dxpy.AppError("Molecules BNX file required for Align Molecules flag (-x)")

        try:
            optargs_xml = dx_utils.download_and_gunzip_file(job_inputs["optargs_xml"])
            scaffold_cmd += ["-q", optargs_xml]

        except KeyError:
            raise dxpy.AppError("OptArgs XML file required for Align Molecules flag (-x)")

    if job_inputs["generate_chimeric"] is True:
        scaffold_cmd += ["-y"]

        if molecules_bnx_file:
            scaffold_cmd += ["-m", molecules_bnx_file]

        else:
            try:
                molecules_bnx_file = dx_utils.download_and_gunzip_file(job_inputs["molecules_bnx_file"])
                scaffold_cmd += ["-m", molecules_bnx_file]

            except KeyError:
                raise dxpy.AppError("Molecules BNX file required for Generate Molecules flag")

        if "err_files" in job_inputs:
            err_files = [dx_utils.download_and_gunzip_file(err_file) for err_file in job_inputs["err_files"]]
            for err in err_files:
                scaffold_cmd += ["-e", err]
    dx_utils.run_cmd(scaffold_cmd)

    dx_utils.run_cmd(["tree", output_dir])

    scaffold_final_ncbi = glob.glob(
        os.path.join(output_dir, 'hybrid_scaffolds', '*_HYBRID_SCAFFOLD_NCBI.fasta'))[0]
    unscaffolded_final = glob.glob(
        os.path.join(output_dir, 'hybrid_scaffolds', '*_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta'))[0]
    scaffold_final = glob.glob(
        os.path.join(output_dir, 'hybrid_scaffolds', '*_HYBRID_SCAFFOLD.fasta'))
    scaffold_final.extend(glob.glob(
        os.path.join(output_dir, 'hybrid_scaffolds', '*_HYBRID_SCAFFOLD.cmap')))
    scaffold_final.extend(glob.glob(
        os.path.join(output_dir, 'hybrid_scaffolds', '*_HYBRID_SCAFFOLD.agp')))

    scaffold_output = glob.glob(os.path.join(output_dir, 'hybrid_scaffolds', '*_HYBRID_SCAFFOLD.xmap'))
    scaffold_output.extend(glob.glob(os.path.join(output_dir, 'hybrid_scaffolds', '*_HYBRID_SCAFFOLD_q.cmap')))
    scaffold_output.extend(glob.glob(os.path.join(output_dir, 'hybrid_scaffolds', '*_HYBRID_SCAFFOLD_r.cmap')))
    scaffold_output = [f for f in scaffold_output if f not in scaffold_final]

    cut_and_conflict = glob.glob(os.path.join(output_dir, 'hybrid_scaffolds*', 'conflicts*.txt'))
    cut_and_conflict.extend(glob.glob(os.path.join(output_dir, 'hybrid_scaffolds*', '*_annotations.bed')))

    # make sure output files don't have colons
    dx_utils.run_cmd(["sed", "-i.bak", "s/:/_/g", scaffold_final_ncbi])
    dx_utils.run_cmd(["sed", "-i.bak", "s/:/_/g", unscaffolded_final])

    # upload outputs
    output = {"scaffold_final": [dx_utils.gzip_and_upload(f) for f in scaffold_final],
            "scaffold_output": [dx_utils.gzip_and_upload(f) for f in scaffold_output],
            "cut_and_conflict": [dxpy.dxlink(dxpy.upload_local_file(f)) for f in cut_and_conflict],
            "ncbi_scaffold_final": dx_utils.gzip_and_upload(scaffold_final_ncbi),
            "unscaffolded_final": dx_utils.gzip_and_upload(unscaffolded_final)}

    tar_name = "hybrid_scaffold_output.tar.gz"
    tar_cmd = "tar czvf {tar_name} {outdir}".format(
        tar_name=tar_name,
        outdir=output_dir)
    dx_utils.run_cmd(tar_cmd)
    output_id = dxpy.upload_local_file(tar_name)
    
    output["scaffold_targz"] = dxpy.dxlink(output_id)
    return output

