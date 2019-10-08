#!/usr/bin/env python
import dxpy, dx_utils, os, subprocess, re, matplotlib, pylab, multiprocessing
matplotlib.use('Agg')
import pylab

GENERATOR_FILENAME = 'generators.txt'
MAX_KMER_COVERAGE = 1000
NUM_LINES_TO_SNIFF = 1000
GENOMESCOPE_MIN_COL = 30
GENOMESCOPE_MAX_COL = 48

def _run_cmd(cmd):
    print subprocess.list2cmdline(cmd)
    subprocess.check_call(cmd)


def plot_hist(hist_fn, output_prefix, mer_length):
    """Create a histogram of the kmer counts"""
    hist_data = {}
    for l in open(hist_fn):
        num_kmers, counts = map(int, l.strip().split())
        hist_data[num_kmers] = counts

    full_hist_data = [hist_data.get(i, 0) for i in xrange(1, num_kmers+1)]

    ofn = '{0}.pdf'.format(output_prefix)
    fig = pylab.figure()
    pylab.bar(xrange(1, len(full_hist_data)+1), full_hist_data, width=1, facecolor='blue', edgecolor='blue')
    pylab.title('Histogram of {0} k-mer counts'.format(mer_length))
    pylab.xlabel('k-mer coverage')
    pylab.ylabel('number of unique k-mers')
    fig.savefig(ofn, dpi=300)

    return ofn


def _get_read_length():
    with open(GENERATOR_FILENAME) as fh:
        cmd = fh.readline().strip()
        cmd += ' | head -n {0}'.format(NUM_LINES_TO_SNIFF)
        fastx_lines = subprocess.check_output(cmd, shell=True)
        if fastx_lines[0] == '>':
            return len(re.sub('\n', '', fastx_lines[1:].partition('>')[0]))
        elif fastx_lines[0] == '@':
            return len(fastx_lines.strip().split('\n')[1])


def _get_sequence_stream(dxf):
    """From the given dxfile, create a command to stream the contents
    to stdout and  gunzip it if needed."""
    fn = dxpy.describe(dxf)['name']
    cmd = 'dx cat {0} '.format(dxf['$dnanexus_link'])
    if os.path.splitext(fn)[-1] == '.gz':
        cmd += '| gunzip '

    return cmd


def _write_generator_file(dxfiles, generator_fn):
    """Creates a file with a line describing how to extract sequences
    from each input dxfile"""
    with open(generator_fn, 'w') as fh:
        for dxf in dxfiles:
            fh.write(_get_sequence_stream(dxf) + '\n')


def _get_genomescope_summary(fn, model_exists):
    genomescope_summary = {}
    if model_exists:
        with open(fn) as fh:
            for l in fh:
                min_val = l[GENOMESCOPE_MIN_COL:GENOMESCOPE_MAX_COL].strip()
                max_val = l[GENOMESCOPE_MAX_COL:].strip()
                if l.startswith('Heterozygosity'):
                    genomescope_summary['genomescope_heterozygosity_estimate'] = '{0} - {1}'.format(min_val, max_val)
                elif l.startswith('Genome Haploid Length'):
                    genomescope_summary['genomescope_haploid_length_estimate'] = '{0} - {1}'.format(min_val, max_val)
                elif l.startswith('Genome Repeat Length'):
                    genomescope_summary['genomescope_repeat_length_estimate'] = '{0} - {1}'.format(min_val, max_val)
                elif l.startswith('Genome Unique Length'):
                    genomescope_summary['genomescope_unique_length_estimate'] = '{0} - {1}'.format(min_val, max_val)
                elif l.startswith('Model Fit'):
                    genomescope_summary['genomescope_model_fit'] = '{0} - {1}'.format(min_val, max_val)
                elif l.startswith('Read Error Rate'):
                    genomescope_summary['genomescope_read_error_rate'] = '{0} - {1}'.format(min_val, max_val)
    else:
        genomescope_summary = {
            'genomescope_heterozygosity_estimate': 'GenomeScope failed to converge',
            'genomescope_haploid_length_estimate': 'GenomeScope failed to converge',
            'genomescope_repeat_length_estimate': 'GenomeScope failed to converge',
            'genomescope_unique_length_estimate': 'GenomeScope failed to converge',
            'genomescope_model_fit': 'GenomeScope failed to converge',
            'genomescope_read_error_rate': 'GenomeScope failed to converge'
        }

    return genomescope_summary

def _run_meryl(output, sequences, k_mer_size, min_k_mer_count):


    meryl_kmer = ["meryl", "count", "threads={}".format(multiprocessing.cpu_count()), "k={}".format(k_mer_size)]
    for file_ref in sequences:

        dx_utils.download_and_gunzip_file(file_ref)
        meryl_kmer.extend([dxpy.describe(file_ref)["name"].replace(".gz",'')])

    meryl_kmer.extend(["output", "out"])
    
    dx_utils.run_cmd(meryl_kmer)
    with open("kmers",'w') as kmers:
        
        subprocess.check_call(["meryl", "print", "greater-than", str(min_k_mer_count), "out"],shell=False,stdout=kmers) 

    with open("mer_counts.tsv",'w') as kmers_index:

        subprocess.check_call(['awk', r'{printf "%s %s\n", NR,$2}', 'kmers'],shell=False, stdout=kmers_index)    

    output["histogram"] = dxpy.dxlink(dxpy.upload_local_file("mer_counts.tsv"))


@dxpy.entry_point('main')
def main(**job_inputs):

    output = {}

    # make k-mer histogra w/ meryl
    _write_generator_file(job_inputs['sequences_fastx'], GENERATOR_FILENAME)
    _run_meryl(output, job_inputs["sequences_fastx"], job_inputs["k_mer_size"], job_inputs["min_k_mer_count"])

    # Run Genomescope
    read_length = _get_read_length()
    cmd = ['Rscript', './genomescope.R', "mer_counts.tsv", str(job_inputs['k_mer_size']),
           str(read_length), './', str(MAX_KMER_COVERAGE)
    ]
    _run_cmd(cmd)
    genomescope_summary = _get_genomescope_summary('summary.txt', os.path.exists('model.txt'))

    # Upload the output files.
    output['genomescope_figures'] = [
        dxpy.upload_local_file('plot.png', name='{0}.gs.png'.format(job_inputs['output_prefix'])),
        dxpy.upload_local_file('plot.log.png', name='{0}.gs.log.png'.format(job_inputs['output_prefix']))
    ]
    output['genomescope_files'] = [
        dxpy.upload_local_file('summary.txt', name='{0}.summary.txt'.format(job_inputs['output_prefix'])),
        dxpy.upload_local_file('progress.txt', name='{0}.progress.txt'.format(job_inputs['output_prefix']))
    ]
    # If GenomeScope failed to converge for some reason, there will be no model.txt file.
    if os.path.exists('model.txt'):
        output['genomescope_files'].append(dxpy.upload_local_file('model.txt', name='{0}.model.txt'.format(job_inputs['output_prefix'])))
    output.update(genomescope_summary)

    return output

dxpy.run()
