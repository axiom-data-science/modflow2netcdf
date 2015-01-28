# coding=utf-8

import os
import tempfile
import textwrap
from datetime import datetime

import ConfigParser

from pyproj import Proj, transform
import numpy as np

from pygc import great_circle

from flopy.modflow.mf import Modflow
import flopy.utils.binaryfile as flopy_binary

import netCDF4
import pytz
from dateutil.parser import parse as date_parse

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
                fills.append(self.lpf.hdry)
            else:
                logger.warning("No LPF file available, using default 'hdry' value of -1e30")
                fills.append(-1e30)
        return fills

    def get_coordinates(self):
        y, x, z = self.dis.get_node_coordinates()

        # Compute grid from top left corner, so switch Y back around
        y = y[::-1]

        # Compute grid from top left corner, so switch X back around
        x = x[::-1]

        # Make Z negative (down)
        z = z*-1.
        # Convert Z to cartesian
        z = z[:, :, ::-1]

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

        # Find output files we want to process.  These are the defaults if
        # no specific file is specified in the configuration
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
        with LoggingTimer("Loading head file", logger.info):
            try:
                if self.head_file:
                    head_obj = flopy_binary.HeadFile(self.head_file, precision=self.precision)
                elif headfile:
                    head_obj = flopy_binary.HeadFile(headfile[1], precision=self.precision)
                else:
                    logger.warning("No Headfile found")
            except BaseException:
                logger.exception("Exception occured when trying to load the HeadFile into flopy. Skipping!")

        # Cell budget file
        cell_obj = None
        with LoggingTimer("Loading cell budget file", logger.info):
            try:
                if self.cbud_file:
                    cell_obj = flopy_binary.CellBudgetFile(self.cbud_file, precision=self.precision)
                elif cellfile:
                    cell_obj = flopy_binary.CellBudgetFile(cellfile[1], precision=self.precision)
                else:
                    logger.warning("No CellBudget file found")
            except BaseException:
                logger.exception("Exception occured when trying to load the CellBudgetFile into flopy. Skipping!")

        return dict(head_obj=head_obj,
                    cell_obj=cell_obj)

    def to_plot(self, variable=None, level=None, time=None, colormap=None):

        import matplotlib.cm as cm
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.axes3d import Axes3D

        if level is not None:
            try:
                level = int(level)
            except ValueError:
                logger.error("Level must be an integer, defaulting to 0.")
                level = 0
        else:
            level = 0

        if time is not None:
            try:
                time = int(time)
            except ValueError:
                logger.error("Time must be an integer, defaulting to 0.")
                time = 0
        else:
            time = 0

        try:
            _, tmp_file = tempfile.mkstemp(suffix=".nc")
            self.to_netcdf(output_file=tmp_file)
            nc = netCDF4.Dataset(tmp_file)

            if variable is not None and variable not in nc.variables:
                raise ValueError("Variable {0} was not found in NetCDF file.  Available variables are: {1}".format(variable, ", ".join([v for v in nc.variables])))

            # Common variables
            x = nc.variables.get("longitude")[:]
            y = nc.variables.get("latitude")[:]
            z = nc.variables.get("elevation")[level, :]

            fig = plt.figure(figsize=(20, 10))

            if colormap is None:
                colormap = cm.Reds

            def plot_thing(rows, columns, spot, camera_height, camera_azimuth):
                ax = fig.add_subplot(rows, columns, spot, projection='3d')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.view_init(camera_height, camera_azimuth)
                if variable is None:
                    ax.set_title('Time: {0} Level: {1}'.format(time, level))
                    p = ax.plot_surface(x, y, z, rstride=1, cstride=1, linewidth=0, cmap=colormap, alpha=0.80)
                    fig.colorbar(p, shrink=0.7)
                else:
                    data = nc.variables.get(variable)
                    data = data[time, level, :, :]
                    m = cm.ScalarMappable(cmap=colormap)
                    m.set_array(data)
                    colors = m.to_rgba(data)[:, :, 0]
                    ax.set_title('Time: {0} Level: {1} Variable: {2}'.format(time, level, variable))
                    p = ax.plot_surface(x, y, z, rstride=1, cstride=1, linewidth=0, alpha=0.80, facecolors=colormap(colors))
                    fig.colorbar(m, shrink=0.7)

            plot_thing(2, 2, 1, 40, 30)
            plot_thing(2, 2, 2, 40, 120)
            plot_thing(2, 2, 3, 40, 210)
            plot_thing(2, 2, 4, 40, 300)

            return plt

        except:
            raise
        finally:
            os.remove(tmp_file)

    def to_netcdf(self, output_file):

        fillvalue = -9999.9

        outputs = self.get_output_objects()
        head_obj = outputs["head_obj"]
        cell_obj = outputs["cell_obj"]

        if os.path.exists(output_file):
            os.remove(output_file)

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
            lat = nc.createVariable('latitude', 'f8', ('x', 'y',), chunksizes=(x_chunk, y_chunk,))
            lat.units         = "degrees_north"
            lat.standard_name = "latitude"
            lat.long_name     = "latitude"
            lat.axis          = "Y"
            lat[:]            = self.ys

            # Longitude
            lon = nc.createVariable('longitude', 'f8', ('x', 'y',), chunksizes=(x_chunk, y_chunk,))
            lon.units         = "degrees_east"
            lon.standard_name = "longitude"
            lon.long_name     = "longitude"
            lon.axis          = "X"
            lon[:]            = self.xs

            # Elevation
            ele = nc.createVariable('elevation', 'f8', ('layer', 'x', 'y',), chunksizes=(z_chunk, x_chunk, y_chunk,), zlib=True)
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

            delc = nc.createVariable('delc', 'f4', ('x',))
            delc.units = 'meters'
            delc.long_name = "Column spacing in the rectangular grid"
            if self.grid_units == 'feet':
                delc[:] = self.dis.delc.array[::-1] * 0.3048
            else:
                delc[:] = self.dis.delc.array[::-1]
            if self.grid_rotation != 0:
                delc.comments = textwrap.dedent("""This is the column spacing that applied to the UNROTATED grid. \
                                This grid HAS been rotated before being saved to NetCDF. \
                                To compute the unrotated grid, use the origin point and this array.""")

            delr = nc.createVariable('delr', 'f4', ('y'))
            delr.units = 'meters'
            delr.long_name = "Row spacing in the rectangular grid"
            if self.grid_units == 'feet':
                delr[:] = self.dis.delr.array[::-1] * 0.3048
            else:
                delr[:] = self.dis.delr.array[::-1]
            if self.grid_rotation != 0:
                delr.comments = textwrap.dedent("""This is the row spacing that applied to the UNROTATED grid. \
                                This grid HAS been rotated before being saved to NetCDF. \
                                To compute the unrotated grid, use the origin point and this array.""")

            # Workaround for CF/CDM.
            # http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/StandardCoordinateTransforms.html
            # "explicit_field"
            exp  = nc.createVariable('VerticalTransform', 'S1')
            exp.transform_name           = "explicit_field"
            exp.existingDataField        = "elevation"
            exp._CoordinateTransformType = "vertical"
            exp._CoordinateAxes          = "layer"

        def create_time(netcdf_file, time_values, t_chunk):
            nc.createDimension("time", len(time_values))
            time = nc.createVariable("time",    "f8", ("time",), chunksizes=(t_chunk,))
            time.units          = "{0} since {1}".format(self.time_units, self.base_date.isoformat().split('.')[0] + "Z")
            time.standard_name  = "time"
            time.long_name      = "time of measurement"
            time.calendar       = "gregorian"
            time[:] = np.asarray(time_values)

        def create_variable(netcdf_file, name, attributes, t_chunk):
            var = nc.createVariable(name, 'f4', ('time', 'layer', 'x', 'y',), fill_value=fillvalue, chunksizes=(t_chunk, z_chunk, x_chunk, y_chunk,), zlib=True)
            var.units         = "{0}^3/time".format(self.grid_units)
            var.standard_name = standard_var_name
            var.long_name     = standard_var_name.upper()
            var.coordinates   = "time layer latitude longitude"
            return var

        # Headfile
        if head_obj is not None:

            with LoggingTimer("Writing HEADS to file", logger.info):
                # Time
                time_values = head_obj.get_times()
                t_chunk = min(len(time_values), 100)
                create_time(nc, time_values, t_chunk)

                head = nc.createVariable('heads', 'f4', ('time', 'layer', 'x', 'y',), fill_value=fillvalue, chunksizes=(t_chunk, z_chunk, x_chunk, y_chunk,), zlib=True)
                head.units         = "{0}^3/time".format(self.grid_units)
                head.standard_name = "heads"
                head.long_name     = "heads"
                head.coordinates   = "time layer latitude longitude"

                for i, time in enumerate(time_values):
                    head_array = head_obj.get_data(totim=time)
                    for f in self.fills:
                        head_array[head_array == f] = fillvalue
                    head_array[head_array <= -1e15] = fillvalue
                    head[i, :, :, :] = head_array[:, :, ::-1]

                head_obj.close()

        if cell_obj is not None:
            if nc.variables.get("time") is None:
                # Time
                time_values = cell_obj.get_times()
                t_chunk = min(len(time_values), 100)
                create_time(nc, time_values, t_chunk)

            for j, var_name in enumerate(cell_obj.unique_record_names()):
                standard_var_name = var_name.strip().lower().replace(' ', '_')
                with LoggingTimer("Writing {} to file".format(standard_var_name), logger.info):
                    attrs = dict(units="{0}^3/time".format(self.grid_units),
                                 standard_name=standard_var_name,
                                 long_name=standard_var_name.upper().replace("_", ""),
                                 coordinates="time layer latitude longitude")
                    var = create_variable(nc, standard_var_name, attrs, t_chunk)

                    # Average the flows onto the grid center
                    centered_variable = None
                    if standard_var_name in ["flow_right_face", "flow_front_face", "flow_lower_face"]:
                        vname = '{0}_centered'.format(standard_var_name)
                        if standard_var_name == "flow_right_face":
                            standard_name = "grid_directed_groundwater_velocity_in_the_u_direction"
                        elif standard_var_name == "flow_front_face":
                            standard_name = "grid_directed_groundwater_velocity_in_the_w_direction"
                        elif standard_var_name == "flow_lower_face":
                            standard_name = "grid_directed_groundwater_velocity_in_the_v_direction"
                        attrs = dict(units="{0}^3/time".format(self.grid_units),
                                     standard_name=standard_name,
                                     long_name=vname.upper().replace("_", ""),
                                     coordinates="time layer latitude longitude")
                        centered_variable = create_variable(nc, vname, attrs, t_chunk)
                        centered_variable[:] = fillvalue

                    for i, time in enumerate(time_values):
                        cell_array = cell_obj.get_data(text=var_name, totim=time, full3D=True)
                        if isinstance(cell_array, list) and len(cell_array) == 1:
                            cell_array = cell_array[0]
                        elif isinstance(cell_array, list) and len(cell_array) > 1:
                            # Sum values if there are more than one in a grid cell.
                            # This was suggested by Chris L. on the 12/22/14 call.
                            cell_array = np.sum(cell_array, axis=0)
                        elif isinstance(cell_array, list) and len(cell_array) == 0:
                            # No data returned
                            logger.warning("No data returned for '{!s}' at time index {!s} ({!s}).".format(var_name.strip(), i, time))
                            cell_array = np.ma.zeros(var.shape[1:])
                            cell_array.mask = True

                        for f in self.fills:
                            cell_array[cell_array == f] = fillvalue
                        cell_array[cell_array <= -1e15] = fillvalue
                        cell_array[cell_array == 0.] = fillvalue
                        var[i, :, :, :] = cell_array[:, :, ::-1]

                        # Average the flows onto the grid center
                        if centered_variable is not None:
                            z, m, n = cell_array.shape
                            new_cell_data = var[i, :, :, :]
                            if standard_var_name == "flow_right_face":
                                # We lose the last column.
                                averaged_array = 0.5 * (new_cell_data[:, :, 0:n-1] + new_cell_data[:, :, 1:n])
                                centered_variable[i, :, :, 0:-1] = averaged_array
                            elif standard_var_name == "flow_lower_face":
                                # We lose the last row
                                averaged_array = 0.5 * (new_cell_data[:, 0:m-1, :] + new_cell_data[:, 1:m, :])
                                centered_variable[i, :, 0:-1, :] = averaged_array
                            elif standard_var_name == "flow_front_face":
                                # We lose the first vertical
                                averaged_array = 0.5 * (new_cell_data[0:z-1, :, :] + new_cell_data[1:z, :, :])
                                centered_variable[i, 1:, :, :] = averaged_array

            cell_obj.close()

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
            self.grid_x        = config.getfloat('space', 'origin_x')
            self.grid_y        = config.getfloat('space', 'origin_y')
            self.grid_rotation = config.getfloat('space', 'rotation')

            # CRS
            config_crs = config.getint('space', 'crs')
            try:
                self.grid_crs = Proj(init='epsg:{0!s}'.format(config_crs))
            except RuntimeError:
                raise ValueError("Could not understand EPSG code '{!s}' from config file.".format(config_crs))

            # Units
            grid_units = config.get('space', 'units')
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

            # Unit of time
            self.time_units = config.get('time', 'units')

            # Base unit of time
            try:
                base_date = date_parse(config.get('time', 'base'))
            except ValueError:
                raise ValueError("Could not parse the base date '{!s}'.".format(base_date))
            else:
                if base_date.tzinfo is None:
                    logger.warning("No timezone information could be extracted from the base date. Assuming UTC.")
                    base_date = base_date.replace(tzinfo=pytz.utc)
                self.base_date = base_date.astimezone(pytz.utc)

        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as e:
            raise ValueError('Configuration file is missing a required section: {0}'.format(e.message))

        # Output files (optional)
        try:
            self.cbud_file = os.path.join(os.path.dirname(config_file), config.get('output', 'cbud'))
            self.head_file = os.path.join(os.path.dirname(config_file), config.get('output', 'head'))
        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as e:
            self.cbud_file = None
            self.head_file = None
