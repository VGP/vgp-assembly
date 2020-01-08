#!/usr/bin/python2.7
import argparse
import os
import re
import sys
import subprocess
import tempfile
import json
from datetime import datetime

import dxpy

BEFORE_FILENAME = 'before-sorted.txt'
AFTER_FILENAME = 'after-sorted.txt'

FILES_TO_IGNORE = {'/var/cache/ldconfig', '/var/lib/dpkg/lock', '/var/lib/dpkg/triggers/Lock'}

###########
# Code from create_asset_trusty
###########

def create_before_file_list():
    before_file_list_path = os.path.join(tempfile.gettempdir(), BEFORE_FILENAME)
    get_system_snapshot(before_file_list_path, [])

def get_file_list(output_file, resources_to_ignore=["/home/dnanexus*"]):
    """
    This method find all the files in the system and writes it to the output file
    """
    tmp_dir = os.path.dirname(output_file) + "*"
    skipped_paths = ["/proc*", tmp_dir, "/run*", "/boot*", "/sys*",
                     "/dev*", "/var/log*", "/root*"]
    cmd = ["sudo", "find", "/"]
    for ignore_dir in skipped_paths + resources_to_ignore:
        cmd.extend(["-not", "-path", ignore_dir])

    env = os.environ.copy()
    env['LC_ALL'] = 'C'
    null_file = open(os.devnull, 'w')
    ps_pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=null_file)
    ps_file = subprocess.Popen(["sort"], stdin=ps_pipe.stdout, stdout=subprocess.PIPE,
                               env=env)

    with open(output_file, "w") as bfile:
        for line in ps_file.stdout:
            sp_code = ps_file.poll()
            file_name = line.rstrip()
            if file_name == "":
                if sp_code is not None:
                    break
                else:
                    continue
            if file_name == "/":
                continue
            try:
                mtime = str(os.path.getmtime(file_name))
            except OSError as os_err:
                #print os_err
                mtime = ''
            # file_name should not have special characters
            bfile.write(file_name + "\t" + str(mtime) + '\n')
    ps_file.stdout.close()

def get_system_snapshot(output_file_path, ignore_files):
    tmp_file_path = tempfile.mktemp()
    get_file_list(tmp_file_path, ignore_files)
    with open(output_file_path, 'w') as output_file_handle:
        proc = subprocess.Popen(['sort', tmp_file_path], stdout=output_file_handle)
        proc.communicate()

def get_file_diffs(first_file, second_file, diff_file):
    """ Get difference between two txt files and write the difference to the
    third file.
    """
    cmd = ["sudo", "comm", "-13", first_file, second_file]
    ps_pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    with open(diff_file, "w") as bfile:
        for line in ps_pipe.stdout:
            line = line.rstrip()
            file_name = '\t'.join(line.split('\t')[:-1])
            if len(file_name) > 0 and file_name not in FILES_TO_IGNORE:
                bfile.write(file_name + '\n')
                print file_name
    ps_pipe.stdout.close()

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('task', choices=['asset', 'snapshot'])
    parser.add_argument(
        "--name",
        required=False,
        default=datetime.now().strftime('%B_%d_%Y_%H_%M'),
        help="Name of the asset or snapshot")
    parser.add_argument(
        "--title",
        required=False,
        help="Title of the asset")
    parser.add_argument(
        "--description",
        required=False,
        help="Description of the asset")
    parser.add_argument(
        "--version",
        required=False,
        help="Version of the asset to build")

    args = parser.parse_args()

    #if args.name is None:
    #    args.name = raw_input('Asset name: ')
    #if args.title is None:
    #    args.title = raw_input('Asset title: ')
    #if args.description is None:
    #    args.description = raw_input('Asset description: ')
    #if args.version is None:
    #    args.version = raw_input('Asset version: ')

    return parser.parse_args()


def get_dist_and_release():
    cmd = 'lsb_release -a '
    distribution = subprocess.check_output(cmd + "| grep Distributor | awk '{print $3}'", shell=True).strip()
    release      = subprocess.check_output(cmd + "| grep Release | awk '{print $2}'", shell=True).strip()

    return (distribution, release)


def create_snapshot(name, ignore_files, suffix='.snapshot', visibility='visible'):
    before_file_list_path = os.path.join(tempfile.gettempdir(), BEFORE_FILENAME)
    after_file_list_path = os.path.join(tempfile.gettempdir(), AFTER_FILENAME)
    get_system_snapshot(after_file_list_path, ignore_files)


    diff_file_path = os.path.join(tempfile.gettempdir(), "diff.txt")
    get_file_diffs(before_file_list_path, after_file_list_path, diff_file_path)

    tar_output = re.sub(r"\s+", '-', name) + suffix
    tar_output = dxpy.PROJECT_CONTEXT_ID + ":" + tar_output
    executable = dxpy.describe(dxpy.JOB_ID)['executable']
    tar_details = {'executable': executable, 
                   'executable_name': dxpy.describe(executable)['name']}
    tar_cmd = ["tar", "-Pcz", "--no-recursion", "-T", diff_file_path, "-f", "-"]
    tar_ps = subprocess.Popen(tar_cmd, stdout=subprocess.PIPE)
    upload_ps = subprocess.Popen(["dx", "upload", "-", "--wait", "--brief", "-o", tar_output, 
                                  "--details", json.dumps(tar_details), "--visibility", visibility],
                                  stdin=tar_ps.stdout, stdout=subprocess.PIPE)

    tar_ps.stdout.close()
    asset_tarball_id = upload_ps.communicate()[0].rstrip()
    tar_ps.wait()
    upload_ps.stdout.close()

    return asset_tarball_id


def create_asset(name, title, description, version):
    asset_tarball_id = create_snapshot(name, ["/home/dnanexus*"], suffix='.tar.gz', visibility='hidden')

    record_name = name
    record_details = {"archiveFileId": {"$dnanexus_link": asset_tarball_id}}
    record_properties = {"version": version,
                         "title": title,
                         "description": description}
    asset_bundle = dxpy.new_dxrecord(name=record_name,
                                     types=["AssetBundle"], details=record_details,
                                     properties=record_properties, close=True,
                                     project=dxpy.PROJECT_CONTEXT_ID)

    with open('dxasset.json', 'w') as fh:
        distribution, release = get_dist_and_release()
        fh.write(json.dumps({'name': name, 'title': title, 'description': description, 'version': version, 'distribution': distribution, 'release': release}))

    return asset_bundle.get_id()

if __name__ == '__main__':
    args = get_parser()
    if args.task == 'asset':
        if args.title is None or args.description is None or args.version is None:
            print 'Assets require a title, description, and version.'
            sys.exit(1)
        create_asset(args.name, args.title, args.description, args.version)
    elif args.task == 'snapshot':
        create_snapshot(args.name, ["/home/dnanexus/dx*"])
