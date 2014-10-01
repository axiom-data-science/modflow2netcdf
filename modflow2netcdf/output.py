# coding=utf-8

import os
import math
from copy import copy
from datetime import datetime

from pyproj import Geod, Proj, transform
import numpy as np
import cartopy

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ioff()

from pygc import great_circle

from flopy.modflow.mf import Modflow

import netCDF4

from modflow2netcdf import logger
from modflow2netcdf.utils import LoggingTimer


class ModflowOutput(object):

    def __init__(self, namfilename, geofile=None, version=None, exe_name=None, verbose=None, model_ws=None):
        if exe_name is None:
            exe_name = 'mf2005.exe'
        if version is None:
            version = 'mf2k'
        if verbose is not True:
            verbose = False

        self.log = get_log(verbose)

        with LoggingTimer("Loaded input and output files", self.log):
            self.mf = Modflow.load(namfilename,
                                   version=version,
                                   exe_name=exe_name,
                                   verbose=verbose,
                                   model_ws=model_ws)
            assert self.mf is not None

        # Load DIS
        self.dis = self.mf.get_package('DIS')
        assert self.dis is not None

        self.get_coordinates(geofile)

    def get_coordinates(self, geofile):
        y, x, z = self.dis.get_node_coordinates()
        #y = y[::-1]

        # Make Z negative (down)
        z = z*-1.
        # Do we need to convert the y axis to cartisian?  I believe we do...
        z = z[:,:,::-1]

        grid_crs, grid_x, grid_y, grid_rotation, grid_units = self.parse_geofile(geofile)
        try:
            provided_crs = Proj(init='epsg:{0!s}'.format(grid_crs))
        except RuntimeError:
            raise ValueError("Could not understand EPSG code '{!s}' from geofile.".format(grid_crs))

        # Convert distances to grid cell centers (from origin) to meters
        if grid_units == 'feet':
            x = x*0.3048
            y = y*0.3048

        # Transform to a known CRS
        wgs84_crs        = Proj(init='epsg:4326')
        known_x, known_y = transform(provided_crs, wgs84_crs, grid_x, grid_y)
        self.log("Input origin point: {!s}, {!s} (EPSG:{!s})".format(grid_x, grid_y, grid_crs))
        self.log("Lat/Lon origin point: {!s}, {!s} (EPSG:4326)".format(known_x, known_y))

        notrotated_xs = np.ndarray(0)
        notrotated_ys = np.ndarray(0)
        upper = great_circle(distance=x, latitude=known_y, longitude=known_x, azimuth=90)
        for top_x, top_y in zip(upper["longitude"], upper["latitude"]):
            # Compute the column points for each point in the upper row.
            # Because this is a regular grid (rectangles), we can just rotate by 180 degrees plus the rotation angle.
            row = great_circle(distance=y, latitude=top_y, longitude=top_x, azimuth=180)
            notrotated_xs = np.append(notrotated_xs, row["longitude"])
            notrotated_ys = np.append(notrotated_ys, row["latitude"])
        # Reshape
        notrotated_xs = notrotated_xs.reshape(self.dis.nrow, self.dis.ncol)
        notrotated_ys = notrotated_ys.reshape(self.dis.nrow, self.dis.ncol)

        if grid_rotation != 0:
            rotated_xs  = np.ndarray(0)
            rotated_ys  = np.ndarray(0)
            upper_rotated   = great_circle(distance=x, latitude=known_y, longitude=known_x, azimuth=90+grid_rotation)
            for top_x, top_y in zip(upper_rotated["longitude"], upper_rotated["latitude"]):
                # Compute the column points for each point in the upper row.
                # Because this is a regular grid (rectangles), we can just rotate by 180 degrees plus the rotation angle.
                row = great_circle(distance=y, latitude=top_y, longitude=top_x, azimuth=180+grid_rotation)
                rotated_xs = np.append(rotated_xs, row["longitude"])
                rotated_ys = np.append(rotated_ys, row["latitude"])
            # Reshape
            rotated_xs = rotated_xs.reshape(self.dis.nrow, self.dis.ncol)
            rotated_ys = rotated_ys.reshape(self.dis.nrow, self.dis.ncol)
        else:
            rotated_ys = notrotated_ys
            rotated_xs = notrotated_xs

        self.origin_x = known_x
        self.origin_y = known_y
        self.no_rotation_xs = notrotated_xs
        self.no_rotation_ys = notrotated_ys
        self.xs = rotated_xs
        self.ys = rotated_ys
        self.zs = z

    def to_plot(self):
        # Setup figure
        fig = plt.figure()
        map_proj = cartopy.crs.PlateCarree()

        vertical_layers = self.zs.shape[0]

        minx = min(np.min(self.xs), np.min(self.no_rotation_xs))-0.025
        miny = min(np.min(self.ys), np.min(self.no_rotation_ys))-0.025
        maxx = max(np.max(self.xs), np.max(self.no_rotation_xs))+0.025
        maxy = max(np.max(self.ys), np.max(self.no_rotation_ys))+0.025

        before = fig.add_subplot(vertical_layers + 1, 2, 1, projection=map_proj)
        before.set_title("Grid")
        before.set_xlim(left=minx, right=maxx)
        before.set_ylim(bottom=miny, top=maxy)
        before.add_feature(cartopy.feature.LAND, edgecolor='black', zorder=0)
        before.pcolor(self.no_rotation_xs, self.no_rotation_ys, np.zeros((self.dis.nrow, self.dis.ncol)))
        before.plot(self.origin_x, self.origin_y, color='red', linewidth=4, marker='x')
        before.text(self.origin_x + 0.005, self.origin_y, 'Origin ({!s}, {!s})'.format(round(self.origin_x, 2), round(self.origin_y, 2)), fontsize=10)

        after = fig.add_subplot(vertical_layers + 1, 2, 2, projection=map_proj)
        after.set_title("Rotated")
        after.set_xlim(left=minx, right=maxx)
        after.set_ylim(bottom=miny, top=maxy)
        after.add_feature(cartopy.feature.LAND, edgecolor='black', zorder=0)
        after.pcolor(self.xs, self.ys, np.zeros((self.dis.nrow, self.dis.ncol)))
        after.plot(self.origin_x, self.origin_y, color='red', linewidth=4, marker='x')
        after.text(self.origin_x + 0.005, self.origin_y, 'Origin ({!s}, {!s})'.format(round(self.origin_x, 2), round(self.origin_y, 2)), fontsize=10)

        # Recalculate the bounds for the images so they reference the rotated grid.
        minx = np.min(self.xs)-0.025
        miny = np.min(self.ys)-0.025
        maxx = np.max(self.xs)+0.025
        maxy = np.max(self.ys)+0.025

        for k in range(vertical_layers):
            after = fig.add_subplot(vertical_layers + 1, 2, k+3, projection=map_proj)
            after.set_title("Vertical Layer {!s}".format(k))
            after.set_xlim(left=minx, right=maxx)
            after.set_ylim(bottom=miny, top=maxy)
            after.add_feature(cartopy.feature.LAND, edgecolor='black', zorder=0)
            after.pcolor(self.xs, self.ys, self.zs[k,:,:])

        plt.show()

    def to_netcdf(self, output_file):
        # Metadata
        t_size = 2
        z_size, x_size, y_size = self.zs.shape
        t_chunk = min(t_size, 100)
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
        nc.createDimension("time", t_size)
        nc.createDimension('x', x_size)
        nc.createDimension('y', y_size)
        nc.createDimension('layer', z_size)

        # Time
        time = nc.createVariable("time",    "f8", ("time",), chunksizes=(t_chunk,))
        time.units          = "seconds since 1970-01-01T00:00:00Z"
        time.standard_name  = "time"
        time.long_name      = "time of measurement"
        time.calendar       = "gregorian"
        time[:] = np.asarray([0, 3600])

        # Metadata variables
        crs = nc.createVariable("crs", "i4")
        crs.long_name           = "http://www.opengis.net/def/crs/EPSG/0/4326"
        crs.epsg_code           = "EPSG:4326"
        crs.semi_major_axis     = float(6378137.0)
        crs.inverse_flattening  = float(298.257223563)

        # Latitude
        lat = nc.createVariable('latitude', 'f8', ('y', 'x',), chunksizes=(y_chunk, x_chunk,))
        lat.units         = "degrees_north"
        lat.standard_name = "latitude"
        lat.long_name     = "latitude"
        lat.axis          = "Y"
        lat[:]            = self.ys

        # Longitude
        lon = nc.createVariable('longitude', 'f8', ('y', 'x',), chunksizes=(y_chunk, x_chunk,))
        lon.units         = "degrees_east"
        lon.standard_name = "longitude"
        lon.long_name     = "longitude"
        lon.axis          = "X"
        lon[:]            = self.xs

        # Elevation
        ele = nc.createVariable('elevation', 'f8', ('layer', 'y', 'x',), chunksizes=(z_chunk, y_chunk, x_chunk,))
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

        # Dummy variable (for now)
        d1 = nc.createVariable('dummy1', 'f4', ('time', 'layer', 'y', 'x',), chunksizes=(t_chunk, z_chunk, y_chunk, x_chunk,))
        d1.units         = "foo"
        d1.standard_name = "bar"
        d1.long_name     = "foo bar"
        d1.coordinates   = "time layer latitude longitude"
        d1[:]            = np.random.random((t_size, z_size, y_size, x_size))

        nc.sync()
        nc.close()
        return nc


    def parse_geofile(self, geofile):
        if geofile is None:
            raise ValueError("No geofile provided")

        try:
            with open(geofile) as f:
                # crs
                # x
                # y
                # rotation
                # units
                info = f.read()
            data = [x.strip() for x in info.split("\n")][0:5]
            data = [x if x != "" and x != "-" else None for x in data]

            try:
                int(data[0])
            except ValueError:
                raise ValueError("CRS code '{!s}' from geofile is not an integer.".format(data[0]))

            if data[4] is None:
                logger.info("Defaulting to 'meters' as the grid units")
                data[4] = "meters"
            elif data[4][0] == "m":
                data[4] = "meters"
            elif data[4][0] == "f":
                data[4] = "feet"
            else:
                raise ValueError("Only units of 'meters' and 'feet' are allowed in the geofile")

            data[3] = float(data[3])

            return data

        except BaseException:
            logger.exception("ok")
            raise ValueError("Could not open geofile: {!s}".format(geofile))


def get_log(verbose=None):
    if verbose:
        return logger.info
    else:
        return logger.debug


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

