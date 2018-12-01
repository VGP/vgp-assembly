"""
This script will upload data to the designated VGP upload bucket in AWS 
and create a link to a corresponding file object in DNAnexus.

The directory structure created in S3 will be as follows:
s3://bucketname/
|--species
|   |-- <species_name>
|   |   |-- <species_id>
|   |   |   |-- genomic_data
|   |   |   |   |-- 10x
|   |   |   |   |-- arima
|   |   |   |   |-- bionano
|   |   |   |   |-- pacbio
|   |   |   |-- transcriptomic_data
|   |   |   |   |-- <tissue_id>
|   |   |   |   |   |-- pacbio
|   |   |   |   |   |-- illumina

In DNAnexus, 'species_id' will be converted to the project name with the directory structure preserved as above

"""

import dxpy
import boto3

import argparse
import os
import subprocess
import sys

import hashlib
from functools import partial

VGP_BUCKET = "genomeark-test"

def parse_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Upload VGP Data')

    ap.add_argument('path',
                    help='filepath or folderpath to upload')

    ap.add_argument('-p', '--profile',
                    help='AWS Profile name',
                    required=True)

    ap.add_argument('-n', '--species-name',
                    help='Species identifier',
                    required=True)

    ap.add_argument('-i', '--species-id',
                    help='Species identifier',
                    required=True)
    ap.add_argument('-t', '--tissue',
                    help='Tissue Type',
                    required=False)
    ap.add_argument('-d', '--datatype',
                    choices=['pacbio', '10x', 'bionano', 'arima', 'illumina', 'phase', 'hic'],
                    help='Sequencing technology or datatype',
                    required=True),
    ap.add_argument('-e', '--local-path',
                    help='file containing local file path, assuming we have the matching aws path. Skip aws s3 sync.',
                    required=False),
    ap.add_argument('-b', '--bill-to',
                    help='Bill to',
                    required=False)

    return ap.parse_args()


def run_cmd(cmd, returnOutput=False):
    print cmd
    if returnOutput:
        output = subprocess.check_output(
            cmd, shell=True, executable='/bin/bash')
        print output
        return output
    else:
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def md5sum(filename):
    with open(filename, mode='rb') as f:
        d = hashlib.md5()
        for buf in iter(partial(f.read, 128), b''):
            d.update(buf)
    return d.hexdigest()

def locate_or_create_dx_drive(drive_name='genomeark'):
    # findDrives is an API method that has not been explicitly added to dxpy yet, so call the API method explicitly.
    drives = dxpy.DXHTTPRequest('/system/findDrives', {'name': drive_name}, always_retry=True)['results']

    # use boto to read the profile info from aws config file
    s3client = boto3.session.Session(profile_name=drive_name)
    profile = s3client.get_credentials()
    
    if len(drives) == 1:
        # Make sure the drive we found is up to date with the latest credentials
        drive_id = drives[0]['id']
        update = {'credentials': {'accessKeyId': profile.access_key,
                                  "secretAccessKey": profile.secret_key}}
        drive_id = dxpy.DXHTTPRequest('/{0}/update'.format(drive_id), update)

        return drive_id
    elif len(drives) == 0:
        # if no drive exists, create it
        print("Creating drive with name: {0}".format(drive_name))

        # create drive using API call
        new_drive_def = {'name': drive_name,
                     'cloud': 'aws',
                     'credentials': {'accessKeyId': profile.access_key,
                                     "secretAccessKey": profile.secret_key}}
        drive_id = dxpy.DXHTTPRequest('/drive/new', new_drive_def)
        print("Created drive with id: {0}".format(drive_id))
        return drive_id
    elif len(drives) > 1:
        print("More than one drives found with name '{0}'".format(drive_name))
        sys.exit(1)


def locate_or_create_dx_project(project_name):
    '''Try to find the project with the given name.  If one doesn't exist,
    we'll create it.'''
    projects = dxpy.find_projects(name=project_name, name_mode='glob', return_handler=True, level='CONTRIBUTE')

    project = [p for p in projects]
    if len(project) < 1:
        project = dxpy.DXProject(dxpy.api.project_new({'name': project_name, 'summary': 'FALCON Unzip Assembly'})['id'])
    elif len(project) > 1:
        print 'Found more than 1 project matching ' + project_name + '.'
        print 'Please provide a unique project!'
        sys.exit(1)
    else:
        project = project[0]

    return project

def main(path, profile, species_name, species_id, datatype, tissue=None, local_path=None, bill_to=None):

    print local_path

    # determine directory to upload to
    if tissue:
        upload_path = os.path.join('species', species_name, species_id, 'transcriptomic_data', tissue, datatype)
    else:
        upload_path = os.path.join('species', species_name, species_id, 'genomic_data', datatype)

    # determine whether file or folder
    if os.path.isfile(path):
        # check size
        # connect to S3 using Boto3 client
        s3client = boto3.session.Session(profile_name=profile).resource('s3')

        # upload file to bucket directly
        object_key = os.path.join(upload_path, os.path.basename(path))

        if local_path:
            # do nothing
	    print 'aws path exists. skip uploading to aws.'
	else:
            s3client.Bucket(VGP_BUCKET).upload_file(path, object_key)

        updated_files = [object_key]
        os_paths = [path]
    elif os.path.isdir(path):
        # here we need to shell out to aws cli to use folder sync
        upload_path = "s3://{0}/{1}".format(VGP_BUCKET, upload_path)

        if local_path:
            # list files in path
            os_paths = open(local_path, 'r').read().splitlines()
	    # Here we don't need to check which files are uploaded. Upload all listed.
            updated_files = os_paths
	    if len(os_paths) > 0:
                os_paths = [f for f in os_paths if f]
                
        else:
            cmd = 'aws s3 sync {0} {1} --no-progress --profile {2}'.format(path, upload_path, profile)
            output = run_cmd(cmd, returnOutput=True)

            # expected output: 
            # upload: test/test.txt to s3://genomeark-test/species/Calypte_anna/bCalAnn1/...
            updated_files = [msg.split(' ')[-1] for msg in output.split('\n')]
            # get rid of empty strings and remove the 's3://bucket-name/' prefix
            updated_files = [f.replace('s3://{0}/'.format(VGP_BUCKET), '') for f in updated_files if f]

            if len(updated_files) > 0:
                os_paths =  [msg.split(' ')[1] for msg in output.split('\n') if len(msg) > 1]
                os_paths = [f for f in os_paths if f]
    else:
        print("Invalid filepath specified")
        sys.exit(1)

    if len(updated_files) == 0:
        print("No new files uploaded")
        sys.exit(0)

    # now link the newly uploaded files to DNAnexus
    dx_drive = locate_or_create_dx_drive(profile)
    dx_project = locate_or_create_dx_project(species_id)

    # in case project properties doesn't have species_name, update it
    dxpy.api.project_set_properties(dx_project.id, input_params={'properties': {'species_name': species_name}})

    for object, os_path in zip(updated_files, os_paths):
        folder_path, filename = os.path.split('/' + object)
        # remove species_name and species_id from folder path
        folder_path = folder_path.replace('/species/{0}/{1}'.format(species_name, species_id), '')
        
        print folder_path, filename
        # create new dx file
        file = dxpy.api.file_new({'project': dx_project.id,
                                'folder': folder_path,
                                'parents': True,
                                'name': filename,
                                'drive': dx_drive["id"],
                                'tags': [species_name, species_id, datatype],
                                'md5sum': md5sum(os_path),
                                'symlinkPath': {
                                    "container": "{region}:{bucket}".format(region='us-east-1', bucket=VGP_BUCKET),
                                    "object": object
                                    },
                                'details': {'container': VGP_BUCKET,
                                            'object': object}
                                  })
        print('Associated object {obj} with DNAnexus file {file}'.format(obj=object, file=file['id']))

if __name__ == '__main__':
    ap = parse_args()
    main(ap.path, ap.profile, ap.species_name, ap.species_id, ap.datatype, ap.tissue, ap.local_path, ap.bill_to)
