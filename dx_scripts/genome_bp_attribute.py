import dxpy
import subprocess
import time
import sys
import ast



def latest_job(name_string,project):
    # --origin-jobs
    job_id = subprocess.check_output('dx find jobs --name {0} --origin-jobs --state done --brief --project {1}'.format(name_string,project),shell=True)
    job_id = job_id.strip().split('\n')
    job_id = filter(None,job_id)
    if len(job_id) == 0:
        return 'job_not_found'
    else:
        return job_id


def fetch_project_id(project_query):
    if len(project_query) == len('project-FQGZx48099bYj52x2PGfJ24q') and project_query.startswith('project-'):
        project_name = subprocess.check_output('dx describe {} --name'.format(project_query),shell=True)
        return project_query,project_name
    else:
        project_id = subprocess.check_output('dx find projects --name {} --json | jq -r .[0].id'.format(project_query), shell=True)
        return project_id,project_query


if __name__ == "__main__":
    project_query=sys.argv[1]
    # get project id
    project_id,project_name=fetch_project_id(project_query)
    genomescope_jobs=latest_job("*Genome*cope*",project_id)
    header = True
    if genomescope_jobs != 'job_not_found':
        for genomescope_job in genomescope_jobs:
            kmer = subprocess.check_output(
                'dx describe {genomescope_job} --json | jq .input.mer_length'.format(genomescope_job=genomescope_job),
                shell=True)
            output = subprocess.check_output(
                'dx describe {genomescope_job} --json | jq .output'.format(genomescope_job=genomescope_job), shell=True)
            output=ast.literal_eval(output)
            if 'genomescope_heterozygosity_estimate' in output:
                if header:
                    print('project_name\tkmer\theterozygosity\trepeat\thaploid')
                    header = False

                heterozygosity,repeat,haploid=output['genomescope_heterozygosity_estimate'],\
                                              output['genomescope_repeat_length_estimate'],\
                                              output['genomescope_haploid_length_estimate']
                print('{project_name}\t{kmer}\t{heterozygosity}\t{repeat}\t{haploid}'.format(project_name=project_name,
                                                                           kmer=kmer,
                                                                           heterozygosity=heterozygosity,
                                                                           repeat=repeat,
                                                                           haploid=haploid))
            else:
                summary_files=subprocess.check_output('dx find data --name summary.txt --brief', shell=True)
                summary_files = summary_files.strip().split('\n')
                summary_files = filter(None, summary_files)
                for summary_file in summary_files:
                    subprocess.check_call('dx cat {}'.format(summary_file), shell=True)
                break

    else:
        print('{project_name\tjobnotfound'.format(project_name=project_name))


