"""
Before running netcdf_cmdline:
    1) Add geographic location information to name file
    2) Change the name of your head file to <namefile_base_name>.hds
    3) Change the extension of your budget file to <namefile_base_name>.cbc

NetcdfCmdline python command line application.  Provides command line access
to Flopy's the NetCDF export functionality.

    Command Line Parameters
    ----------
    -n : File Path
        Path to MODFLOW name file cooresponding to the project being exported
        to NetCDF.
    -ni : File Path
        Path to NetCDF output file that will store model input data
    -no : File Path
        Path to NetCDF output file that will store model output data
    -p : String ('single' or 'double')
        Precision of modflow output files, single or double.

"""

import os
import argparse
import flopy

class NetCDFCmdline():
    """
    NetCDFCmdline utility class.  Reads the command line for this script.

    Attributes
    ----------
    name_file_path : string
        Path to MODFLOW name file cooresponding to the project being exported
        to NetCDF.
    model_in_netcdf : string
        Path to NetCDF output file that will store model input data
    model_out_netcdf : string
        Path to NetCDF output file that will store model output data
    precision : string
        Precision of modflow output files, single or double.
    """
    def __init__(self):
        # parse command line
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('-n', '--name_file',
                            help="Modflow Namefile (.nam) to load.",
                            required=True)
        self.parser.add_argument('-ni', '--netcdf_modelinput_file',
                            help="Name of output netcdf file that stores model input.",
                            required=True)
        self.parser.add_argument('-no', '--netcdf_modeloutput_file',
                            help="Name of output netcdf file that stores model output",
                            default='output.nc',
                            required=True)
        self.parser.add_argument('-p', '--precision',
                            help="Precision of modflow output files (single or double)",
                            default='single',
                            required=True)
        self.args = self.parser.parse_args()

        # format and store command line variables
        self.name_file_path   = self.args.name_file
        self.model_in_netcdf = os.path.realpath(self.args.netcdf_modelinput_file)
        self.model_out_netcdf = os.path.realpath(self.args.netcdf_modeloutput_file)
        self.precision = self.args.precision


def run():
    """
    Command line application's main entry point.
    """
    # load command line arguments
    cmdline = NetCDFCmdline()

    # change working directory to model directory
    original_dir = os.getcwd()

    try:
        os.chdir(os.path.dirname(cmdline.name_file_path))
    except:
        print 'Invalid path to name file.'

    try:
        name_file_fullname = os.path.basename(cmdline.name_file_path)

        # load modflow model
        ml = flopy.modflow.Modflow.load(name_file_fullname, check=False)

        # export model inputs
        ml.export(cmdline.model_in_netcdf)

        # export outputs
        flopy.export.utils.output_helper(cmdline.model_out_netcdf,ml,ml.load_results(as_dict=True,precision=cmdline.precision))
    finally:
        # clean up
        os.chdir(original_dir)


if __name__ == "__main__":
    run()

