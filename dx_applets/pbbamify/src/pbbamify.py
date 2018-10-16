import os
import glob
import dxpy
import dx_utils
import tempfile
import multiprocessing

def _get_filenames(files):
    input = {'objects': [file['$dnanexus_link'] for file in files]}
    descriptions = dxpy.DXHTTPRequest('/system/describeDataObjects', input)

    fns = [d['describe']['name'] for d in descriptions['results']]

    return fns

def _make_fofn(filenames):
    with tempfile.NamedTemporaryFile(delete=False, suffix='.fofn', dir='/home/dnanexus') as fh:
        fh.write('\n'.join(filenames))
    
    return fh.name

@dxpy.entry_point("main")
def main(**job_inputs):
    # Download the inputs we can't stream in
    subread_bam_files = [dx_utils.download_and_gunzip_file(f) for f in job_inputs['subreads_bams']]
    ref_fa = dx_utils.download_and_gunzip_file(job_inputs['reference_fa'])

    # create fofn of bam files
    bam_fofn = _make_fofn(subread_bam_files)
    help_cmd = ('pbbamify -h')
    dx_utils.run_cmd(help_cmd)
    
    # pbindex the bam files
    for bam_fn in subread_bam_files:
        pbindex_cmd = 'pbindex {0}'.format(bam_fn)
        dx_utils.run_cmd(pbindex_cmd)
    
    # stream the input and output
    output_pbbams = []
    output_pbbais = []
    for in_dxlink in job_inputs['aligned_bams']:
        bam_fn = dx_utils.download_and_gunzip_file(in_dxlink)
        out_bam_fn = bam_fn.replace('.bam', '.tmp.bam')
        sorted_fn = bam_fn.replace('.bam', '.pb.bam')

        pbbamify_cmd = 'pbbamify --input {bam_fn} --output {bam_out} {ref_fa} {bam_fofn} '
        pbbamify_cmd = pbbamify_cmd.format(bam_fn = bam_fn, bam_out = out_bam_fn,
                                           ref_fa=ref_fa, bam_fofn=bam_fofn)
        dx_utils.run_cmd(pbbamify_cmd)
        
        cmd = 'samtools sort -m 2G -@ {0} "{1}" -o "{2}" '.format(
        multiprocessing.cpu_count(), out_bam_fn, sorted_fn)
        dx_utils.run_cmd(cmd)

        cmd = 'samtools index {0}'.format(sorted_fn)
        dx_utils.run_cmd(cmd)

        output_pbbams.append(dxpy.dxlink(dxpy.upload_local_file(sorted_fn)))
        output_pbbais.append(dxpy.dxlink(dxpy.upload_local_file(sorted_fn + '.bai')))
    return {"pb_bam": output_pbbams}
