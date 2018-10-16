#!/usr/bin/env python
from __future__ import print_function
import unittest
import yaml
import json
import subprocess
import dxpy

TEST_DATA_YAML = './test_data.yml'
MAIN_APP = '../'


def remove_output(output_value):
    if type(output_value) is list:
        map(remove_output, output_value)
    elif dxpy.is_dxlink(output_value):
        dxpy.remove(output_value)
    elif type(output_value) is dict:
        map(remove_output, output_value.values())


class Minimap2Tests(unittest.TestCase):
    def setUp(self):
        print('Building applets')
        with open(TEST_DATA_YAML) as fh:
            self.test_data = self._filter_test_data(yaml.safe_load(fh))

        self.main_applet = json.loads(subprocess.check_output(['dx', 'build', '-f', MAIN_APP]).strip())['id']


    def tearDown(self):
        print('Removing applets')
        subprocess.check_call(['dx', 'rm', self.main_applet])
        print('Removing job outputs.')
        map(lambda job: remove_output(job.describe()['output']), self.jobs)


    @staticmethod
    def _filter_test_data(test_data):
        region = dxpy.describe(dxpy.PROJECT_CONTEXT_ID)['region']
        env = None
        if dxpy.APISERVER_HOST.startswith('staging'):
            test_data = test_data['staging'].get(region, None)
            env = 'staging'
        else:
            test_data = test_data['production'].get(region, None)
            env = 'production'

        if test_data is None:
            raise Exception('Could not find test description for {0} in the {1} region.'.format(env, region))

        return test_data


    @staticmethod
    def get_sequencing_technology_choices():
        for i in json.load(open('../dxapp.json'))['inputSpec']:
            if i['name'] == 'sequencing_technology':
                return i['choices']

        raise Exception('Could not sequencing_technology input argument in dxapp.json.')


    def test_applet(self):
        self.assertIsNotNone(self.test_data)

        main_applet = dxpy.DXApplet(self.main_applet)
        self.jobs = []
        sequencing_technologies = self.get_sequencing_technology_choices()
        for sequencing_technology in sequencing_technologies:
            job_input = {'reads_fastqgz': dxpy.dxlink(self.test_data['reads_fastqgz']),
                         'reads2_fastqgz': dxpy.dxlink(self.test_data['reads2_fastqgz']),
                         'genome_fastagz': dxpy.dxlink(self.test_data['genome_fastagz']),
                         'sequencing_technology': sequencing_technology}
            job = main_applet.run(job_input)
            print('Running {0} test in {1}'.format(sequencing_technology, job))
            self.jobs.append(job)

        for job in self.jobs:
            job.wait_on_done()
