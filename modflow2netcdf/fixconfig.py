#-------------------------------------------------------------------------------
# Name:        FixConfig
# Purpose:     Fix configuration file to work with mfnetcdf.  This includes
#              transforming origin_x and origin_y into the crs used by mfnetcdf.
#
# Author:      spaulins
#
# Created:     25/03/2015
# Copyright:   (c) spaulins 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------


import ConfigParser
import argparse
import os
from pyproj import Proj, transform  # cartopgraphic transformations

def main():
    # Parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--config_file',
                        help="Modflow2netcdf configiration file to load.",
                        required=True)
    parser.add_argument('-r', '--crs_code',
                        help="Modflow2netcdf configiration file to load.",
                        required=True)

    args = parser.parse_args()
    config_file_path = os.path.realpath(args.config_file)
    crs_code = args.crs_code
    config_file = FixConfigFile(config_file_path)
    config_file.fixCRS(crs_code)
    config_file.update_config_file()

class FixConfigFile(object):
    def __init__(self, config_file_path):
        self.config_file_path = config_file_path

        ## Parse config file ##
        if self.config_file_path is None:
            raise ValueError("No config file provided")

        self.config = ConfigParser.SafeConfigParser()
        try:
            self.config.read(config_file_path)
        except ConfigParser.ParsingError as e:
            raise ValueError("Bad configuration file.  Please check the contents. {!s}.".format(e.message))

    def update_config_file(self):
        # Write updated config file
        with open(self.config_file_path, 'wb') as configfile:
            self.config.write(configfile)

    def fixCRS(self, crs_code):
        # CRS
        config_crs = self.config.getint('space', 'crs')
        origin_x   = self.config.getfloat('space', 'origin_x')
        origin_y   = self.config.getfloat('space', 'origin_y')
        if config_crs != crs_code:
            # Transform origin_x and origin_y to correct projection
            try:
                grid_crs = Proj(init='epsg:{0!s}'.format(config_crs))
            except RuntimeError:
                raise ValueError("Could not understand EPSG code '{!s}' from config file.".format(self.config_crs))

            wgs84_crs        = Proj(init='epsg:%s' % crs_code)
            transform_x, transform_y = transform(grid_crs, wgs84_crs, origin_x, origin_y)
            print transform_x

            self.config.set('space', 'crs', crs_code)
            self.config.set('space', 'origin_x', str(transform_x))
            self.config.set('space', 'origin_y', str(transform_y))

if __name__ == "__main__":
    main()
