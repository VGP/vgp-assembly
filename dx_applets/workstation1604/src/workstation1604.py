#!/usr/bin/env python
# cloud_workstation 1.0.0

import os
import re
import subprocess
import dxpy

def _run_cmd(cmd):
    print (cmd)
    subprocess.check_call(cmd, shell=True)

@dxpy.entry_point('main')
def main(**job_inputs):
    if 'fids' in job_inputs:
        for fid in job_inputs['fids']:
            cmd = 'dx download {0}'.format(fid['$dnanexus_link'])
            _run_cmd(cmd)

    seconds_to_sleep = dxpy.utils.normalize_timedelta(job_inputs['max_session_length'])/1000
    cmd = 'sleep {0}'.format(seconds_to_sleep)
    _run_cmd(cmd)
    output = {}

    return output

dxpy.run()