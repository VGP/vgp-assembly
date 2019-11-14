from __future__ import print_function
import dxpy
import subprocess
import time


app_id_hardcode_up_version={
    "app-FPkkQ4j0gjx97J1X1496B9zF":"purge_haplotig (bitbucket v1.0.3+ 1.Nov.2018)",
    "applet-FZzYjYQ0j3b9pVP06qg4Q97Y": "purge_dups github ca23030ccf4254dfd2d3a5ea90d0eed41c24f88b",
    "applet-Fb721zj0j3b1xV916kB8Y14J": "purge_dups github ca23030ccf4254dfd2d3a5ea90d0eed41c24f88b",
    "app-FVVgJ7Q09zJpb0KZ9P3v1BpQ": "Solve3.2.1_04122018",
    "app-FVpb0j00px8fVZ9qPPGYxxP8": "Salsa 2.2",
    "app-FXF87GQ0yV32v3Q32v06xBvv": "smrtlink_6.0.0.47841",
    "app-Fb0JBK8012x8z3gG91Yxyj3q": "smrtlink_7.0.1.66975",
    "app-FPgQ4Y8086pf03z5J04ZkXF3": "longranger 2.2.2",
    "applet-FZY5j400j3bP4b62GxKB057v": "freebayes 1.3.1",
    "applet-FZB2yG80j3b75K260z2yF688": "freebayes 1.3.1",
    "applet-FYVY3X00x9Y0YpKZ6xgjQxZG": "freebayes v1.2.0-17-ga78ffc0 parallel by chr",
}

def latest_job(name_string):

    if '||' not in name_string:
        job_id = subprocess.check_output('dx find jobs --name {0} --all-jobs --state done -n 1 --brief'.format(name_string),shell=True)
        return job_id.strip()
    else:
        name_list = name_string.split('||')
        for name_list_member in name_list:
            job_id = subprocess.check_output('dx find jobs --name {0} --all-jobs --state done -n 1 --brief'.format(name_list_member), shell=True)
            if job_id != '':
                return job_id.strip()


def job_2_app(job_id):
    try:
        app_id = dxpy.describe(job_id)['app']

    except KeyError:
        app_id = dxpy.describe(job_id)['applet']
    return app_id.strip()

def app_2_version(app_id):
    try:
        version = dxpy.describe(app_id)['version']
    except KeyError:
        version = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(dxpy.describe(app_id)['created']))
    return version

def app_2_upversion(app_id):
    try:
        upversion = dxpy.describe(app_id)['details']['upstreamVersion']
    except KeyError:
        if app_id in app_id_hardcode_up_version:
            return app_id_hardcode_up_version[app_id]
        else:
            upversion = 'NA'
    return upversion.strip()

def start_time(job_id):
    epoch_time=dxpy.describe(job_id)['startedRunning']
    startedRunning = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(epoch_time))
    return startedRunning

falcon_job_id=latest_job("*alcon*aligner*")
falcon_unzip_job_id=latest_job("*nzip*olish*||*Unzip*Haplotype*Arrow")
purge_job_id=latest_job("*purge*")
scaff10x_job_id=latest_job("Scaff10x*")
bionano_job_id=latest_job("Bionano*")
salsa_job_id=latest_job("Salsa*")
polish_job_id=latest_job("polish")
longranger_job_id=latest_job("10X*Longranger*Align*")
freebayes_job_id=latest_job("Free*ayes*")



print('\t'.join(['falcon','falcon_unzip','purge','scaff10x','bionano','salsa','polish','longranger','freebayes']))

job_list = [falcon_job_id,falcon_unzip_job_id,purge_job_id,scaff10x_job_id,bionano_job_id,
              salsa_job_id,polish_job_id,longranger_job_id,freebayes_job_id]
print('\t'.join(job_list))

app_list = map(job_2_app,job_list)
print('\t'.join(map(start_time,job_list)))
print('\t'.join(app_list))
print('\t'.join(map(app_2_version,app_list)))
print('\t'.join(map(app_2_upversion,app_list)))

