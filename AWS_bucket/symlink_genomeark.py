import dxpy
import boto3

import argparse
import os
import subprocess
import sys

import hashlib
from functools import partial

VGP_BUCKET = "genomeark"

def parse_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Symlink VGP data to dnanexus')

    ap.add_argument('object',
                    help='s3 object to symlink from genomeark (s3 path)')
    
    ap.add_argument('-l', '--local',
                    help='local file path to generate md5sum',
                    required=True)

    ap.add_argument('-p', '--profile',
                    help='AWS Profile name',
                    required=True)

    ap.add_argument('-n', '--species-name',
                    help='Species name',
                    required=True)

    ap.add_argument('-i', '--species-id',
                    help='Species identifier',
                    required=True)

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

def main(object, local, profile, species_name, species_id):

    # set profile and project
    dx_drive = locate_or_create_dx_drive(profile)
    dx_project = locate_or_create_dx_project(species_id)

    # in case project properties doesn't have species_name, update it
    dxpy.api.project_set_properties(dx_project.id, input_params={'properties': {'species_name': species_name}})

    # format object if includes full path
    obj = object.replace('s3://genomeark/species/{0}/{1}'.format(species_name, species_id), '')
    obj = obj.strip('/')
    obj = os.path.join('species', species_name, species_id, obj)

    # remove species_name and species_id from folder path to match dnanexus
    folder_path = obj.replace('species/{0}/{1}'.format(species_name, species_id), '')
    folder_path = folder_path.replace('/{0}'.format(local),'')
    filename = os.path.basename(obj)
    
    # Debugging #        
    print ('Path on DNAnexus: {0}/{1}'.format(folder_path, filename))

    # create new dx file
    file = dxpy.api.file_new({'project': dx_project.id,
                            'folder': folder_path,
                            'parents': True,
                            'name': filename,
                            'drive': dx_drive["id"],
                            'tags': [species_name, species_id],
                            'md5sum': md5sum(local),
                            'symlinkPath': {
                                "container": "{region}:{bucket}".format(region='us-east-1', bucket=VGP_BUCKET),
                                "object": obj
                                },
                            'details': {'container': VGP_BUCKET,
                                        'object':  obj}
                              })
    print('Associated object {obj} with DNAnexus file {file}'.format(obj=obj, file=file['id']))

if __name__ == '__main__':
    ap = parse_args()
    # object, local, profile, species_name, species_id
    main(ap.object, ap.local, ap.profile, ap.species_name, ap.species_id)

