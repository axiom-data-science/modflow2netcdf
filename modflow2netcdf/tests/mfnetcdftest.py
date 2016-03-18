#-------------------------------------------------------------------------------
# Name:        TestOutput
# Purpose:     Perform unit tests on ModflowOutput
#
# Author:      Scott Paulinski
#
# Created:     25/03/2015
#-------------------------------------------------------------------------------

import os
import unittest
from data_validation_engine import ValidationLog
from unit_test_data_util import ProjectUnitTestData, NetCDFTestProject, ExternalNetCDFTestProject
from mfnetcdf import ModflowToNetCDF, VerifyException


class TestModflowToNetCDF(unittest.TestCase):
    # Success of all unit tests
    success = True

    @classmethod
    def setUpClass(self):
        print 'setup'
        # Constants
        self.test_project_list_file = 'TestProjectList.txt'
        self.test_project_log_file = 'TestProjectLog.txt'

        self.log_file = ValidationLog(self.test_project_log_file)

        # Load list of test projects from config file (unittestdatautil)
        self.internal_test_projects = []
        self.external_test_projects = []
        with open(self.test_project_list_file) as project_list_file:
            project_list_file.next()
            for line in project_list_file:
                arrline = line.split(',')
                if arrline[1]:
                    self.external_test_projects.append(ExternalNetCDFTestProject(arrline[0], arrline[2], arrline[3],
                                                                                 arrline[4], arrline[5], arrline[6],
                                                                                 arrline[7], int(arrline[8]),
                                                                                 int(arrline[9]), int(arrline[10])))
                else:
                    self.internal_test_projects.append(NetCDFTestProject(arrline[0], arrline[2], arrline[3], arrline[4],
                                                                         arrline[5]))
            project_list_file.close()

    @classmethod
    def tearDownClass(self):
        self.log_file.write_success(TestModflowToNetCDF.success)
        self.log_file.close()

    def test_internal(self):
        print 'internal'
        return self._test_loop(self.internal_test_projects)

    def test_external(self):
        print 'external'
        return self._test_loop(self.external_test_projects)

    def _test_loop(self, test_projects):
        local_success = True

        for test_project in test_projects:
            try:
                # Run ModflowToNetCDF in verify mode for current project
                fdtest = open(test_project.name_file)
                fdtest.close()
                mfnetcdf = ModflowToNetCDF(test_project.name_file, config_file=test_project.mf2netcdf_config,
                                           exe_name="mf-2005", verbose=True, model_ws=test_project.project_path)
                mfnetcdf.save_netcdf(test_project.mf2netcdf_output_file, verify=True)
            except VerifyException as ex:
                # Verification error occurred
                self.log_file.write_test_failure(test_project.project_name, ex.message)
                local_success = False
            except Exception as ex:
                # An exception other than a verification error was raised during execution
                self.log_file.write_test_failure(test_project.project_name, ex.message)
                local_success = False

            print 'Running tests for %s' % (test_project.name_file)
            if not test_project.run_tests(self.log_file):
                local_success = False

        if not local_success:
            TestModflowToNetCDF.success = False
            raise unittest.TestCase.failureException()

        return local_success

if __name__ == '__main__':
    unittest.main()
