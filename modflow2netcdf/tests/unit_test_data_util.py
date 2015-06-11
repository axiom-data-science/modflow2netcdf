#-------------------------------------------------------------------------------
# Name:        UnitTestDataUtil
# Purpose:
#
# Author:      spaulins
#
# Created:     25/03/2015
# Copyright:   (c) spaulins 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import data_validation_engine
import os
import netCDF4

class NetCDFTestProject(object):
    def __init__(self, project_name, project_path, name_file, mf2netcdf_config, mf2netcdf_output_file):
        self.current_dir = os.getcwd()

        self.project_name = project_name
        self.name_file = os.path.join(project_path, name_file)
        self.mf2netcdf_config =  os.path.join(project_path, mf2netcdf_config)
        self.project_path = project_path
        self.mf2netcdf_output_file = mf2netcdf_output_file

    def run_tests(self, log_file):
        return True

# Configuration information for a specific test project
class ExternalNetCDFTestProject(NetCDFTestProject):
    def __init__(self, project_name, project_path, name_file, mf2netcdf_config, mf2netcdf_output_file,
                 data_verification_config_file, data_verification_file_path, layers, rows, columns):
        # Call base class constructor
        super(ExternalNetCDFTestProject, self).__init__(project_name, project_path, name_file, mf2netcdf_config,
                                                        mf2netcdf_output_file)

        self.data_verification_config_file = os.path.join(self.current_dir, data_verification_config_file)
        self.data_verification_file_path = os.path.join(self.current_dir, data_verification_file_path)
        self.layers = layers
        self.rows = rows
        self.columns = columns

    def run_tests(self, log_file):
        print 'running external tests'
        project_data = ProjectUnitTestData(self.data_verification_config_file, self.data_verification_file_path, log_file)

        # Read from netcdf file
        ncdf = netCDF4.Dataset(self.mf2netcdf_output_file)

        # Loop through data names in the netcdf file
        for varname in ncdf.variables:
            print 'testing %s' % varname
            # Get data
            data = ncdf.variables.get(varname)

            # Verify data
            project_data.confirm_netcdf_result(varname, data)

        project_data.record_missing_data()


class TestProjectVerificationData(object):
    def __init__(self, name, data_check, data_type, data_file_path, data_file_name, data_delimiter, data_error_threshold):
        self.name = name
        self.data_check = data_check
        self.data_type = data_type
        self.data_file_path = data_file_path
        self.data_file_name = data_file_name
        self.data_delimiter = data_delimiter
        self.data_error_threshold = float(data_error_threshold)


# Purpose:  Tests your program's output data against unit test data for a specific project.
# Input:    objProject - TestDataProject
class ProjectUnitTestData(object):
    def __init__(self, config_file, data_files_path, log_file):
        self.project_data = {}
        self.data_test_complete = {}
        self.log_file = log_file

        # Optional Logging
        if self.log_file is not None:
            self.log_file.write_test_header(config_file)

        # Read project config file
        with open(config_file, 'r') as data_verification_config_file:
            data_verification_config_file.next()
            for line in data_verification_config_file:
                arrline = line.split(',')
                complete_path = os.path.join(data_files_path, arrline[3])
                self.project_data[arrline[0]] = TestProjectVerificationData(arrline[0], arrline[1], arrline[2], complete_path, arrline[4], arrline[5], arrline[6])
            data_verification_config_file.close()

    def confirm_netcdf_result(self, result_name, results):
        # Figure out if result fits into any data set defined in the config file based on the result name
        if result_name in self.project_data:
            # Optional Logging
            if self.log_file is not None:
                self.log_file.write_test_data_header(result_name)

            # Check NetCDF data vs. expected results
            tester = data_validation_engine.test_engine_factory(self.project_data[result_name].data_type,
                                                                self.project_data[result_name].data_check,
                                                                self.project_data[result_name].data_file_path,
                                                                self.log_file,
                                                                self.project_data[result_name].data_delimiter,
                                                                self.project_data[result_name].data_error_threshold)
            if tester is None:
                # Report that data is not being tested due to no functional support for that data type
                if self.log_file is not None:
                    self.log_file.write_data_type_not_supported(self.project_data[result_name].data_type)
                    return -1

            self.data_test_complete[result_name] = tester.validate(self.project_data[result_name].data_file_name, results, self.project_data[result_name].data_delimiter)
            if self.log_file is not None:
                self.log_file.write_test_result(self.data_test_complete[result_name], result_name)
        else:
            # Record as data not in expected results
            self.log_file.write_test_data_without_test(result_name)

    def record_missing_data(self):
        for result_name in self.project_data:
            if  not self.data_test_complete.has_key(result_name):
                self.log_file.write_test_data_missing(result_name)