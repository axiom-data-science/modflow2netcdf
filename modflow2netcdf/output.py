# coding=utf-8

import os
import math
from copy import copy
from datetime import datetime

import ConfigParser

from pyproj import Geod, Proj, transform
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ioff()

from pygc import great_circle

from flopy.modflow.mf import Modflow
import flopy.utils.binaryfile as flopy_binary

import netCDF4

from modflow2netcdf import logger
from modflow2netcdf.utils import LoggingTimer


class ModflowOutput(object):

    def __init__(self, namfilename, config_file=None, version=None, exe_name=None, verbose=None, model_ws=None):
        if exe_name is None:
            exe_name = 'mf2005.exe'
        if version is None:
            version = 'mf2k'
        if verbose is not True:
            verbose = False

        with LoggingTimer("Loaded input and output files", logger.info):
            self.mf = Modflow.load(namfilename,
                                   version=version,
                                   exe_name=exe_name,
                                   verbose=verbose,
                                   model_ws=model_ws)
            assert self.mf is not None

        # DIS (required)
        try:
            self.dis = self.mf.get_package('DIS')
            assert self.dis is not None
        except AssertionError:
            raise ValueError("DIS file could not be found.  It is required.")

        # BAS
        try:
            self.bas = self.mf.get_package('BAS6')
            assert self.bas is not None
        except AssertionError:
            logger.warning("BAS file could not be found.")

        # LPF
        try:
            self.lpf = self.mf.get_package('LPF')
            assert self.bas is not None
        except AssertionError:
            logger.warning("LPF file could not be found.")

        # Now process
        with LoggingTimer("Parsing config file", logger.info):
            self.parse_config_file(config_file)
        self.get_coordinates()

        self.fills = self.get_fill_values()

    def get_fill_values(self):

        fills = []
        # Load first fill value from .bas file
        with LoggingTimer("Retrieving fill values", logger.info):
            # hnoflo
            if self.bas is not None:
                fills.append(self.bas.hnoflo)
            else:
                logger.warning("No BAS file available, using default 'hnoflow' value of -999.99")
                fills.append(-999.99)
            # hdry
            if self.lpf is not None:
                logger.warning("No LPF file available, using default 'hdry' value of -1e30")
                fills.append(self.lpf.hdry)
            else:
                fills.append(-1e30)
        return fills

    def get_coordinates(self):
        y, x, z = self.dis.get_node_coordinates()

        # Compute grid from top left corner, so switch Y back around
        y = y[::-1]

        # Make Z negative (down)
        z = z*-1.
        # Do we need to convert the z axis to cartisian?  I believe we do...
        z = z[:,:,::-1]

        # Convert distances to grid cell centers (from origin) to meters
        if self.grid_units == 'feet':
            x = x*0.3048
            y = y*0.3048

        # Transform to a known CRS
        wgs84_crs        = Proj(init='epsg:4326')
        known_x, known_y = transform(self.grid_crs, wgs84_crs, self.grid_x, self.grid_y)
        logger.debug("Input origin point: {!s}, {!s}".format(self.grid_x, self.grid_y))
        logger.debug("Lat/Lon origin point: {!s}, {!s} (EPSG:4326)".format(known_x, known_y))

        notrotated_xs = np.ndarray(0)
        notrotated_ys = np.ndarray(0)
        with LoggingTimer("Computing unrotated output grid", logger.info):
            upper = great_circle(distance=x, latitude=known_y, longitude=known_x, azimuth=90)
            for top_x, top_y in zip(upper["longitude"], upper["latitude"]):
                # Compute the column points for each point in the upper row.
                # Because this is a regular grid (rectangles), we can just rotate by 180 degrees plus the rotation angle.
                row = great_circle(distance=y, latitude=top_y, longitude=top_x, azimuth=180)
                notrotated_xs = np.append(notrotated_xs, row["longitude"])
                notrotated_ys = np.append(notrotated_ys, row["latitude"])

        if self.grid_rotation != 0:
            with LoggingTimer("Computing rotated output grid", logger.info):
                rotated_xs  = np.ndarray(0)
                rotated_ys  = np.ndarray(0)
                upper_rotated   = great_circle(distance=x, latitude=known_y, longitude=known_x, azimuth=90+self.grid_rotation)
                for top_x, top_y in zip(upper_rotated["longitude"], upper_rotated["latitude"]):
                    # Compute the column points for each point in the upper row.
                    # Because this is a regular grid (rectangles), we can just rotate by 180 degrees plus the rotation angle.
                    row = great_circle(distance=y, latitude=top_y, longitude=top_x, azimuth=180+self.grid_rotation)
                    rotated_xs = np.append(rotated_xs, row["longitude"])
                    rotated_ys = np.append(rotated_ys, row["latitude"])
        else:
            rotated_ys = notrotated_ys
            rotated_xs = notrotated_xs

        self.origin_x = known_x
        self.origin_y = known_y
        self.no_rotation_xs = notrotated_xs.reshape(self.dis.ncol, self.dis.nrow).T
        self.no_rotation_ys = notrotated_ys.reshape(self.dis.ncol, self.dis.nrow).T
        self.xs = rotated_xs.reshape(self.dis.ncol, self.dis.nrow).T
        self.ys = rotated_ys.reshape(self.dis.ncol, self.dis.nrow).T
        self.zs = z

    def get_output_objects(self):
        headfile = tuple()
        cellfile = tuple()

        # Find output files we want to process
        for u, f, b in zip(self.mf.external_units, self.mf.external_fnames, self.mf.external_binflag):
            _, ext = os.path.splitext(f)
            if ext == ".bud":
                cellfile = (u, f, b)
            elif ext == ".bhd":
                headfile = (u, f, b)
            else:
                pass

        # Headfile
        head_obj = None
        if headfile:
            with LoggingTimer("Loading head file (.bud)", logger.info):
                head_obj = flopy_binary.HeadFile(headfile[1], precision=self.precision)
        else:
            logger.warning("No Headfile found")

        # Cell budget file
        cell_obj = None
        if cellfile:
            with LoggingTimer("Loading cell budget file (.bhd)", logger.info):
                cell_obj = flopy_binary.CellBudgetFile(cellfile[1], precision=self.precision)
        else:
            logger.warning("No Budfile found")

        return dict(head_obj=head_obj,
                    cell_obj=cell_obj)

    def to_plot(self):

        fillvalue = -9999.9

        outputs = self.get_output_objects()
        head_obj = outputs["head_obj"]
        cell_obj = outputs["cell_obj"]

        # Setup figure
        with LoggingTimer("Plotting", logger.info):
            fig = plt.figure()

            vertical_layers = self.zs.shape[0]
            time_slices = 0
            if head_obj:
                times = head_obj.get_times()
                # Show just one timestep per layer (the last)
                time_slices = vertical_layers

            minx = min(np.min(self.xs), np.min(self.no_rotation_xs))-0.025
            miny = min(np.min(self.ys), np.min(self.no_rotation_ys))-0.025
            maxx = max(np.max(self.xs), np.max(self.no_rotation_xs))+0.025
            maxy = max(np.max(self.ys), np.max(self.no_rotation_ys))+0.025

            before = fig.add_subplot(time_slices + 2, 2, 1)
            before.set_title("Grid")
            before.set_xlim(left=minx, right=maxx)
            before.set_ylim(bottom=miny, top=maxy)
            before.pcolor(self.no_rotation_xs, self.no_rotation_ys, np.zeros((self.dis.nrow, self.dis.ncol)))
            before.plot(self.origin_x, self.origin_y, color='red', linewidth=4, marker='x')
            before.text(self.origin_x + 0.005, self.origin_y, 'Origin ({!s}, {!s})'.format(round(self.origin_x, 2), round(self.origin_y, 2)), fontsize=10)

            after = fig.add_subplot(time_slices + 2, 2, 2)
            after.set_title("Rotated")
            after.set_xlim(left=minx, right=maxx)
            after.set_ylim(bottom=miny, top=maxy)
            after.pcolor(self.xs, self.ys, np.zeros((self.dis.nrow, self.dis.ncol)))
            after.plot(self.origin_x, self.origin_y, color='red', linewidth=4, marker='x')
            after.text(self.origin_x + 0.005, self.origin_y, 'Origin ({!s}, {!s})'.format(round(self.origin_x, 2), round(self.origin_y, 2)), fontsize=10)

            # Recalculate the bounds for the images so they reference the rotated grid.
            minx = np.min(self.xs)-0.025
            miny = np.min(self.ys)-0.025
            maxx = np.max(self.xs)+0.025
            maxy = np.max(self.ys)+0.025

            if head_obj:
                for k in range(time_slices):
                    head = fig.add_subplot(time_slices + 2, 2, k+3)
                    head.set_title("Vertical Layer {!s}".format(k))
                    head.set_xlim(left=minx, right=maxx)
                    head.set_ylim(bottom=miny, top=maxy)
                    head_array = head_obj.get_data(totim=times[-1])
                    for f in self.fills:
                        head_array[head_array==f] = fillvalue
                    head_array[head_array<=-1e15] = fillvalue
                    head.pcolor(self.xs, self.ys, head_array[k, :, :])

        plt.show()

    def to_netcdf(self, output_file):

        fillvalue = -9999.9

        outputs = self.get_output_objects()
        head_obj = outputs["head_obj"]
        cell_obj = outputs["cell_obj"]


        with LoggingTimer("Setting up NetCDF file", logger.info):

            # Metadata
            z_size, x_size, y_size = self.zs.shape
            x_chunk = x_size / 2.
            y_chunk = y_size / 2.
            z_chunk = 1
            min_vertical = np.min(self.zs)
            max_vertical = np.max(self.zs)

            nc = netCDF4.Dataset(output_file, "w")
            nc.setncattr("Conventions", "CF-1.6")
            nc.setncattr("date_created", datetime.utcnow().strftime("%Y-%m-%dT%H:%M:00Z"))
            nc.setncattr("geospatial_vertical_positive",   "up")
            nc.setncattr("geospatial_vertical_min",        min_vertical)
            nc.setncattr("geospatial_vertical_max",        max_vertical)
            nc.setncattr("geospatial_vertical_resolution", "variable")
            nc.setncattr("featureType", "Grid")

            # Dimensions
            nc.createDimension('x', x_size)
            nc.createDimension('y', y_size)
            nc.createDimension('layer', z_size)

            # Metadata variables
            crs = nc.createVariable("crs", "i4")
            crs.long_name           = "http://www.opengis.net/def/crs/EPSG/0/4326"
            crs.epsg_code           = "EPSG:4326"
            crs.semi_major_axis     = float(6378137.0)
            crs.inverse_flattening  = float(298.257223563)

            # Latitude
            lat = nc.createVariable('latitude', 'f8', ('x', 'y',), chunksizes=(y_chunk, x_chunk,))
            lat.units         = "degrees_north"
            lat.standard_name = "latitude"
            lat.long_name     = "latitude"
            lat.axis          = "Y"
            lat[:]            = self.ys

            # Longitude
            lon = nc.createVariable('longitude', 'f8', ('x', 'y',), chunksizes=(y_chunk, x_chunk,))
            lon.units         = "degrees_east"
            lon.standard_name = "longitude"
            lon.long_name     = "longitude"
            lon.axis          = "X"
            lon[:]            = self.xs

            # Elevation
            ele = nc.createVariable('elevation', 'f8', ('layer', 'x', 'y',), chunksizes=(z_chunk, x_chunk, y_chunk,))
            ele.units         = "meters"
            ele.standard_name = "elevation"
            ele.long_name     = "elevation"
            ele.valid_min     = min_vertical
            ele.valid_max     = max_vertical
            ele.positive      = "down"
            ele[:]            = self.zs

            lay = nc.createVariable('layer', 'f4', ('layer',))
            lay.units         = ''
            lay.long_name     = 'layer'
            lay.positive      = "down"
            lay.axis          = 'Z'
            lay[:]            = np.arange(0, z_size)

            # Workaround for CF/CDM.
            # http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/StandardCoordinateTransforms.html
            # "explicit_field"
            exp  = nc.createVariable('VerticalTransform', 'S1')
            exp.transform_name           = "explicit_field"
            exp.existingDataField        = "elevation"
            exp._CoordinateTransformType = "vertical"
            exp._CoordinateAxes          = "layer"

            # Headfile
            if head_obj is not None:

                # Time
                times = head_obj.get_times()
                t_size = len(times)
                t_chunk = min(t_size, 100)

                nc.createDimension("time", t_size)
                time = nc.createVariable("time",    "f8", ("time",), chunksizes=(t_chunk,))
                time.units          = "seconds since 1970-01-01T00:00:00Z"
                time.standard_name  = "time"
                time.long_name      = "time of measurement"
                time.calendar       = "gregorian"
                time[:] = np.asarray(times)

                head = nc.createVariable('heads', 'f4', ('time', 'layer', 'x', 'y',), fill_value=fillvalue, chunksizes=(t_chunk, z_chunk, x_chunk, y_chunk,))
                head.units         = "units of head"
                head.standard_name = "heads standard name"
                head.long_name     = "heads long name"
                head.coordinates   = "time layer latitude longitude"

                for i, time in enumerate(times):
                    head_array = head_obj.get_data(totim=time)
                    for f in self.fills:
                        head_array[head_array==f] = fillvalue
                    head_array[head_array<=-1e15] = fillvalue
                    head[i, :, :, :] = head_array

        with LoggingTimer("Writing NetCDF file", logger.info):
            nc.sync()
            nc.close()

        return nc


    def parse_config_file(self, config_file):
        if config_file is None:
            raise ValueError("No config file provided")

        config = ConfigParser.SafeConfigParser()
        try:
            config.read(config_file)
        except ConfigParser.ParsingError as e:
            raise ValueError("Bad configuration file.  Please check the contents. {!s}.".format(e.message))

        try:
            self.grid_x        = config.getfloat('grid', 'origin_x')
            self.grid_y        = config.getfloat('grid', 'origin_y')
            self.grid_rotation = config.getfloat('grid', 'rotation')

            # CRS
            config_crs = config.getint('grid', 'crs')
            try:
                self.grid_crs = Proj(init='epsg:{0!s}'.format(config_crs))
            except RuntimeError:
                raise ValueError("Could not understand EPSG code '{!s}' from config file.".format(config_crs))

            # Units
            grid_units    = config.get('grid', 'units')
            if grid_units is None:
                logger.info("Defaulting to 'meters' as the grid units")
                grid_units = "meters"
            elif grid_units[0] == "m":
                grid_units = "meters"
            elif grid_units[0] == "f":
                grid_units = "feet"
            else:
                raise ValueError("Only units of 'meters' or 'feet' are allowed in the config file")
            self.grid_units = grid_units

            self.precision = config.get('general', 'precision')

        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as e:
            raise ValueError(e.message)


def parse(namfilename, packages):
    '''
    function to parse the nam file and return a dictionary with types, names, units and handles
    '''
    # add the .nam extension to namfilename if missing
    if os.path.splitext(namfilename)[-1] != ".nam":
        namfilename += '.nam'

    # inititate the ext_unit_dict dictionary
    ext_unit_dict = dict()

    logger.info("Parsing the namefile --> {0:s}".format(namefilename))
    logger.info("Setting filehandles:")

    with open(namfilename, 'r') as f:
        for line in f.readlines():
            line = line.strip().split()

            # Skip empty lines
            if not line:
                continue
            # Skip comments
            if line[0] == '#':
                continue
            # Skip lines that don't have an integer as the second value
            try:
                int(line[1])
            except ValueError:
                continue

            # Figure out which mode to open the file in (ASCII or Binary)
            filemode = 'r'
            if line[0].upper() == 'DATA(BINARY)':
                filemode = 'rb'

            filename = os.path.join(os.path.dirname(namfilename), line[2])
            try:
                with open(filename, filemode) as f:
                    # Make sure it is a valid file
                    f.read()
                ext_unit_dict[int(line[1])] = NameData(line[0], filename, filemode, packages)
            except BaseException:
                logger.error("Could not open referenced file --> {0:s}".format(line[2]))
                continue

    return ext_unit_dict

