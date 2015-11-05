# coding=utf-8

# Standard library imports
import math
import os
import tempfile
import textwrap
import ConfigParser
from datetime import datetime
# Inter-Project imports
from modflow2netcdf import logger
from modflow2netcdf.utils import LoggingTimer
# USGS imports
from flopy.modflow.mf import Modflow
import flopy.utils.binaryfile as flopy_binary
import flopy.utils.formattedfile as flopy_formatted
# Third part package imports
import numpy as np
import netCDF4
from pygc import great_circle
from pyproj import Proj, transform  # cartopgraphic transformations


class VerifyException(Exception):
    def __init__(self, message):
        self.message = message

class ModflowToNetCDF(object):
    """
    ModflowNetCDF Class.

    Parameters
    ----------
    namfilename : MODFLOW project name file (full or relative path)
    config_file : ModflowToNetCDF config file (full or relative path)
    version : Version of MODFLOW
    exe_name : Modflow executable name
    verbose : Run ModflowtoNetCDF in verbose mode
    model_ws : Path to model input files (full or relative path)

    Attributes
    ----------

    Methods
    -------
    to_plot :  Creates a plot.
    save_netcdf : Creates a netcdf file.

    See Also
    --------

    Notes
    -----

    Examples
    --------

    >>> mf = ModflowNetCDF('freyberg.nam', config_file='freyberg.geo', exe_name="mf2005", verbose=False, model_ws='.\\tests\\resources\\freyberg')
    >>> freyberg_output = 'temp_freyberg_output.nc'
    >>> mf.to_netcdf(output_file=freyberg_output)

    """
    def __init__(self, namfilename, config_file=None, version=None, exe_name=None, verbose=None, model_ws=None):
        # Constants
        self.c_max_error = 1e-10
        self.c_epsg_code = '4326'

        if exe_name is None:
            exe_name = 'mf2005.exe'
        if version is None:
            version = 'mf2k'
        if verbose is not True:
            verbose = False

        self.netcdfdata = {}
        with LoggingTimer("Loaded input and output files", logger.info):
            try:
                self.mf = Modflow.load(namfilename,
                                       version=version,
                                       exe_name=exe_name,
                                       verbose=verbose,
                                       model_ws=model_ws)
                assert self.mf is not None
            except AssertionError:
                raise ValueError("Modflow project could not be loaded.")

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

        # Process config file
        with LoggingTimer("Parsing config file", logger.info):
                # Try parsing again
                if self._parse_config_file(config_file) == False:
                    raise Exception('Configuration file error.  Unable to transform origin_x and origin_y to crs %s'
                                    'Either transformation is not available or the Pyproj dependency may be '
                                    'missing.' % self.c_epsg_code)

        self._get_coordinates()

        # Set hnoflo and hdry values
        self.fills = self._get_fill_values()
        self.fills.append(-9999)  # This was needed for the carolina model
        logger.info("Fill values set to: {!s}".format(self.fills))

    def to_plot(self, variable=None, level=None, time=None, colormap=None):
        """
        Plot modflow project map

        Parameters
        ----------
        variable : string
            variable name of data to plot
        level : int
            layer to plot at
        time : int
            time to plot data at
        colormap : matplotlib.cm
            colormap to display

        Returns
        ----------
        matplotlib.pyplot
            object containing plotted figures

        """
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
            tmpfd, tmp_file = tempfile.mkstemp(suffix=".nc")
            os.close(tmpfd)
            self.save_netcdf(output_file=tmp_file)
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

            nc.close()
            return plt

        except:
            raise
        finally:
            os.remove(tmp_file)

    def save_netcdf(self, output_file, verify=False):
        # Create a NetCDF time variable
        def _create_time(time_values, t_chunk):
            nc.createDimension("time", len(time_values))
            time = nc.createVariable("time",    "f8", ("time",), chunksizes=(t_chunk,))
            time.units          = "{0} since {1}".format(self.time_units, self.base_date)
            time.standard_name  = "time"
            time.axis = "T"
            time._CoordinateAxisType = "Time"
            time.ancillary_variables = ""
            time.comment = ""
            time.ioos_category = "Time"
            time.long_name      = "time of measurement"
            time.calendar       = "gregorian"
            self.netcdfdata["time"] = np.asarray(time_values)
            time[:] = self.netcdfdata["time"]

        """
        Save modflow project data to netcdf file

        Parameters
        ----------
        output_file : string
            File name including full file path of output netcdf file
        verify : boolean
            Netcdf file verification mode.  When set to True data in the netcdf output file is verified.
        Returns
        ----------
        netcdf4.dataset
            dataset containing the contents of the outputted netcdf file

        """
        self.fillvalue = -9999.9
        self.cbc_fillvalue = 0.0

        if self.grid_units == 'feet':
            self.netcdfdata["delc"] = self.dis.delc.array[::-1] * 0.3048
            self.netcdfdata["delr"] = self.dis.delr.array[::-1] * 0.3048
        else:
            self.netcdfdata["delc"] = self.dis.delc.array[::-1]
            self.netcdfdata["delr"] = self.dis.delr.array[::-1]

        outputs = self._get_output_objects()
        head_obj = outputs["head_obj"]
        cell_obj = outputs["cell_obj"]

        if os.path.exists(output_file):
            try:
                os.remove(output_file)
            except:
                logger.exception("Could not remove existing NetCDF output file %s.".format(output_file))
                raise Exception("Could not remove existing NetCDF output file %s.".format(output_file))

        with LoggingTimer("Setting up NetCDF file", logger.info):
            # Metadata
            z_size, x_size, y_size = self.netcdfdata["elevation"].shape
            self.netcdfdata["layer"] = np.arange(0, z_size)

            x_chunk = int(x_size / 4) + 1
            y_chunk = int(y_size / 4) + 1
            z_chunk = z_size
            min_vertical = np.min(self.netcdfdata["elevation"])
            max_vertical = np.max(self.netcdfdata["elevation"])

            # Build NetCDF file
            nc = netCDF4.Dataset(output_file, "w")
            nc.setncattr("Conventions", "CF-1.6")
            nc.setncattr("date_created", datetime.utcnow().strftime("%Y-%m-%dT%H:%M:00Z"))
            nc.setncattr("geospatial_vertical_positive",   "up")
            nc.setncattr("geospatial_vertical_min",        min_vertical)
            nc.setncattr("geospatial_vertical_max",        max_vertical)
            nc.setncattr("geospatial_vertical_resolution", "variable")
            nc.setncattr("featureType", "Grid")
            nc.setncattr("origin_x", self.grid_x)
            nc.setncattr("origin_y", self.grid_y)
            nc.setncattr("origin_crs", self.c_epsg_code)
            nc.setncattr("grid_rotation_from_origin", self.grid_rotation)
            for k, v in self.global_attributes.iteritems():
                try:
                    nc.setncattr(k, v)
                except:
                    logger.exception("Could not set global attribute {}.  Check that its value is supported in NetCDF4.".format(k))

            # Dimensions
            nc.createDimension('x', x_size)
            nc.createDimension('y', y_size)
            nc.createDimension('layer', z_size)

            # Metadata variables
            crs = nc.createVariable("crs", "i4")
            crs.long_name           = "http://www.opengis.net/def/crs/EPSG/0/%s" % self.c_epsg_code
            crs.epsg_code           = "EPSG:%s" % self.c_epsg_code
            crs.semi_major_axis     = float(6378137.0)
            crs.inverse_flattening  = float(298.257223563)

            # Latitude
            lat = nc.createVariable('latitude', 'f8', ('x', 'y',), chunksizes=(x_chunk, y_chunk,))
            lat.units         = "degrees_north"
            lat.standard_name = "latitude"
            lat.long_name     = "latitude"
            lat.axis          = "Y"
            lat._CoordinateAxisType = "Lat"
            lat[:]            = self.netcdfdata["latitude"]

            # Longitude
            lon = nc.createVariable('longitude', 'f8', ('x', 'y',), chunksizes=(x_chunk, y_chunk,))
            lon.units         = "degrees_east"
            lon.standard_name = "longitude"
            lon.long_name     = "longitude"
            lon.axis          = "X"
            lon._CoordinateAxisType = "Lon"
            lon[:]            = self.netcdfdata["longitude"]

            # Time
            time_values_head = None
            time_values_cell = None
            time_values_all = None
            # Get times from head and cell file
            if head_obj is not None:
                time_values_head = head_obj.get_times()
            if cell_obj is not None:
                time_values_cell = cell_obj.get_times()
            # Merge times from head and cell file
            if time_values_head and time_values_cell:
                time_values_all = time_values_head + list(set(time_values_cell) - set(time_values_head))
            elif time_values_head:
                time_values_all = time_values_head
            elif time_values_cell:
                time_values_all = time_values_cell
            # Build time NetCDF object
            if time_values_all:
                t_chunk = min(len(time_values_all), 100)
                _create_time(time_values_all, t_chunk)

            # Elevation
            ele = nc.createVariable('elevation', 'f8', ('layer', 'x', 'y',), chunksizes=(z_chunk, x_chunk, y_chunk,), zlib=True)
            ele.units         = "meters"
            ele.standard_name = "elevation"
            ele.long_name     = "elevation"
            ele.valid_min     = min_vertical
            ele.valid_max     = max_vertical
            ele.positive      = "down"
            ele[:]            = self.netcdfdata["elevation"]

            # Layers
            lay = nc.createVariable('layer', 'f4', ('layer',))
            lay.units         = ''
            lay.long_name     = 'layer'
            lay.positive      = "down"
            lay.axis          = 'Z'
            lay[:]            = self.netcdfdata["layer"]

            # Column lengths
            delc = nc.createVariable('delc', 'f4', ('x',))
            delc.units = 'meters'
            delc.long_name = "Column spacing in the rectangular grid"
            delc.setncattr("origin_x", self.grid_x)
            delc.setncattr("origin_y", self.grid_y)
            delc.setncattr("origin_crs", self.c_epsg_code)
            delc[:] = self.netcdfdata["delc"]

            if self.grid_rotation != 0:
                delc.comments = textwrap.dedent("""This is the column spacing that applied to the UNROTATED grid. \
                                This grid HAS been rotated before being saved to NetCDF. \
                                To compute the unrotated grid, use the origin point and this array.""")

            # Row lengths
            delr = nc.createVariable('delr', 'f4', ('y'))
            delr.units = 'meters'
            delr.long_name = "Row spacing in the rectangular grid"
            delr.setncattr("origin_x", self.grid_x)
            delr.setncattr("origin_y", self.grid_y)
            delr.setncattr("origin_crs", self.c_epsg_code)
            delr[:] = self.netcdfdata["delr"]

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

        # Create a NetCDF variable
        def _create_variable(name, attributes, t_chunk, single_layer=False, use_fill_value=True):
            # Normalize variable name
            if not single_layer:
                if use_fill_value:
                    var = nc.createVariable(name, 'f4', ('time', 'layer', 'x', 'y',), fill_value=self.fillvalue, chunksizes=(t_chunk, z_chunk, x_chunk, y_chunk,), zlib=True)
                else:
                    var = nc.createVariable(name, 'f4', ('time', 'layer', 'x', 'y',), chunksizes=(t_chunk, z_chunk, x_chunk, y_chunk,), zlib=True)
            else:
                if use_fill_value:
                    var = nc.createVariable(name, 'f3', ('time', 'x', 'y',), fill_value=self.fillvalue, chunksizes=(t_chunk, x_chunk, y_chunk,), zlib=True)
                else:
                    var = nc.createVariable(name, 'f3', ('time', 'x', 'y',), chunksizes=(t_chunk, x_chunk, y_chunk,), zlib=True)
            for k, v in attributes.iteritems():
                try:
                    var.setncattr(k, v)
                except:
                    logger.exception("Could not set variable attribute {} on {}.  Check that its value is supported in NetCDF4.".format(k, name))
            return var

        # Store contents of head file
        if head_obj is not None:
            with LoggingTimer("Writing HEAD to file", logger.info):
                # Time
                time_values = head_obj.get_times()

                attrs = dict(standard_name='head',
                             long_name='head',
                             coordinates='time lat lon alt',
                             source=' ',
                             references=' ',
                             grid_mapping='crs',
                             scale_factor=1.,
                             add_offset=0.,
                             cf_role='timeseries_id',
                             units="{0}".format(self.grid_units))
                head = _create_variable('head', attrs, t_chunk)
                self.netcdfdata["head"] = 0

                # Loop through all each time of available head data
                for i, time in enumerate(time_values):
                    # Read the head data from the head file
                    head_array = head_obj.get_data(totim=time)
                    # Set up fill values
                    for f in self.fills:
                        head_array[head_array == f] = self.fillvalue
                    head_array[head_array <= -1e15] = self.fillvalue
                    # Store in NetCDF variable
                    head[i, :, :, :] = head_array[:, :, ::-1]

        # Store contents of cell by cell flow file
        if cell_obj is not None:
            # Time
            time_values = cell_obj.get_times()

            for j, var_name in enumerate(cell_obj.unique_record_names()):
                standard_var_name = var_name.strip().lower().replace(' ', '_')
                with LoggingTimer("Writing {} to file".format(standard_var_name), logger.info):
                    # Create NetCDF variable to store current record
                    attrs = dict(units="{0}^3/time".format(self.grid_units),
                                 standard_name=standard_var_name,
                                 long_name=standard_var_name.upper().replace("_", ""),
                                 coordinates="time layer latitude longitude")
                    name = standard_var_name.replace('.', '_').replace(' ', '_').replace('-', '_')
                    print 'Processing %s' % (name)
                    var = np.zeros(shape=(len(time_values), z_size, x_size, y_size), dtype=np.float64)
                    var[:] = self.cbc_fillvalue
                    var_dimension = 3
                    self.netcdfdata[name] = ('CELL', var_name, time_values, False)

                    # Create NetCDF variable to store flows averaged onto the grid center
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
                        centered_name = vname.replace('.', '_').replace(' ', '_').replace('-', '_')
                        centered_variable = _create_variable(centered_name, attrs, t_chunk, use_fill_value=False)
                        self.netcdfdata[centered_name] = ('CELL', var_name, time_values, True)
                        centered_variable[:] = self.cbc_fillvalue

                    # Loop through all times in the current record
                    for i, time in enumerate(time_values):
                        cell_array = self._get_cell_array(cell_obj, var_name, time, var.shape[1:])

                        # Store cell data in NetCDF
                        try:
                            if len(cell_array.shape) == 2:
                                # Returned a 2D array, so set data to the top layer
                                var[i, 0, :, :] = cell_array[:, ::-1]
                            else:
                                # Returned a 3D array
                                var_dimension = 4
                                var[i, :, :, :] = cell_array[:, :, ::-1]

                        except IndexError:
                            logger.exception("Could not save variable {!s} into NetCDF file.  Trying to fit {!s} into {!s}".format(var_name.strip(), var[i, :, :, :].shape, cell_array.shape))

                        # Average the flows onto the grid center
                        if centered_variable is not None:
                            new_cell_data = var[i, :, :, :]
                            centered_data = centered_variable[i, :, :, :]
                            self._average_flows(new_cell_data, centered_data, standard_var_name)
                            centered_variable[i, :, :, :] = centered_data[:, :, :]

                    # Write to standard NetCDF variable
                    if var_dimension == 3:
                        var_ncdf = _create_variable(name, attrs, t_chunk, True, use_fill_value=False)
                        var_ncdf[:, :, :] = var[:, 0, :, :]
                        print 'Saved %s as 3-D variable.' % (name)
                    else:
                        var_ncdf = _create_variable(name, attrs, t_chunk, use_fill_value=False)
                        var_ncdf[:, :, :, :] = var[:, :, :, :]
                        print 'Saved %s as 4-D variable.' % (name)

        with LoggingTimer("Writing NetCDF file", logger.info):
            nc.sync()
            nc.close()

        # Read NetCDF file and verify that data in equals data out
        if verify:
            # Test code
            if head_obj is not None:
                head_obj.close()
            if cell_obj is not None:
                cell_obj.close()
            outputs = self._get_output_objects()
            head_obj = outputs["head_obj"]
            cell_obj = outputs["cell_obj"]


            # Verify data in NetCDF file by loading data and comparing it to ModFLOW data
            objNCData = netCDF4.Dataset(output_file)

            # Loop through data names in the netcdf file
            self.data_found = {}
            for strVarName in objNCData.variables:
                #  if data is on verification list
                if self.netcdfdata.has_key(strVarName):
                    # if head data
                    if strVarName == 'head':
                        # Pass data to head data comparison function
                        self._verify_head_data(objNCData, head_obj, self.fillvalue)
                    # else if cell data
                    elif self.netcdfdata[strVarName][0] == 'CELL':
                        # Pass data and var name to cell data comparison function
                        self._verify_cell_data(objNCData, cell_obj, strVarName, self.netcdfdata[strVarName][1],
                                               self.netcdfdata[strVarName][2], self.netcdfdata[strVarName][3])
                    else:
                        # Verify data and mark data as found
                        data_out = objNCData.variables[strVarName][:]
                        if self._array_comp(data_out, self.netcdfdata[strVarName]) == False:
                            raise VerifyException('Data verification error: %s does not match.  Unmatched arrays dumped to debug_array_*.txt' % strVarName)

                    # Mark data as found
                    self.data_found[strVarName] = 1

            # Verify that all data is accounted for
            for name in self.netcdfdata:
                if self.data_found.has_key(name) == False:
                    raise VerifyException('Data verification error.  Expected data missing from NetCDF file: %s' % name)

        if head_obj is not None:
            head_obj.close()
        if cell_obj is not None:
            cell_obj.close()
        return nc

# save_cdf support methods #
    def _get_output_objects(self):
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
                    if self.head_type == 'formatted':
                        head_obj = flopy_formatted.FormattedHeadFile(self.head_file, precision=self.precision)
                    elif self.head_type == 'binary':
                        head_obj = flopy_binary.HeadFile(self.head_file, precision=self.precision)
                    else:
                        logger.warning("Head file type " + self.head_type + " not supported.  Only 'binary'" +
                                       " and 'formatted' types supported.")
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

    def _get_cell_array(self, cell_obj, var_name, time, empty_shape):
        # Get cell data from flopy
        cell_array = cell_obj.get_data(text=var_name, totim=time, full3D=True, verbose=True)
        if isinstance(cell_array, list) and len(cell_array) == 1:
            cell_array = cell_array[0]
        elif isinstance(cell_array, list) and len(cell_array) > 1:
            # Sum values if there are more than one in a grid cell.
            # This was suggested by Chris L. on the 12/22/14 call.
            cell_array = np.sum(cell_array, axis=0)
        elif isinstance(cell_array, list) and len(cell_array) == 0:
            # No data returned
            logger.warning("No data returned for '{!s}' at time {!s}.".format(var_name.strip(), time))
            cell_array = np.zeros(empty_shape)

        return cell_array

    def _average_flows(self, new_cell_data, centered_data, standard_var_name):
        z, m, n = new_cell_data.shape
        if standard_var_name == "flow_right_face":
            # We lose the last column.
            averaged_array = 0.5 * (new_cell_data[:, :, 0:n-1] + new_cell_data[:, :, 1:n])
            centered_data[:, :, 0:-1] = averaged_array
        elif standard_var_name == "flow_lower_face":
            # We lose the last row
            averaged_array = 0.5 * (new_cell_data[:, 0:m-1, :] + new_cell_data[:, 1:m, :])
            centered_data[:, 0:-1, :] = averaged_array
        elif standard_var_name == "flow_front_face":
            # We lose the first vertical
            averaged_array = 0.5 * (new_cell_data[0:z-1, :, :] + new_cell_data[1:z, :, :])
            centered_data[1:, :, :] = averaged_array

    def _array_comp(self, first_array, second_array):
        diff = first_array - second_array
        max = np.max(np.abs(diff))
        if max > self.c_max_error:
            return False
        return True

# save_cdf internal test methods #
    def _print_array_diff(self, first_array, second_array, first_array_name, second_array_name):
        try:
            diff = first_array - second_array
            self._save_array(first_array_name, first_array)
            self._save_array(second_array_name, second_array)
            self._save_array('debug_array_diff.txt', diff)
        except:
            logger.exception("An error occurred while outputting array differences.")
            return False
        return True

    # Saves an array with up to three dimensions
    def _save_array(self, filename, multi_array):
        with file(filename, 'w') as outfile:
            outfile.write('%s\n' % str(multi_array.shape))
            if len(multi_array.shape) == 4:
                for slice in multi_array:
                    for second_slice in slice:
                        np.savetxt(outfile, second_slice, fmt='%10.6e')
                        outfile.write('\n')
                    outfile.write('\n')
            elif len(multi_array.shape) == 3:
                for slice in multi_array:
                    np.savetxt(outfile, slice, fmt='%10.6e')
                    outfile.write('\n')
            else:
                np.savetxt(outfile, multi_array, fmt='%10.6e')

    def _verify_cell_data(self, objNCData, cell_obj, cdf_var_name, flopy_var_name, time_values, averaged):

        # Get cell data from NetCDF
        netcdf_cell = objNCData.variables[cdf_var_name][:]

        # Loop through time values that should be stored in NetCDF
        for i, time in enumerate(time_values):
            # Get cell data from NetCDF file
            if len(netcdf_cell.shape) == 3:
                netcdf_data = netcdf_cell[i, :, :]
            else:
                netcdf_data = netcdf_cell[i, :, :, :]
            # Get cell data for proper variable name and current time value from flopy
            cell_array = self._get_cell_array(cell_obj, flopy_var_name, time, netcdf_cell.shape[1:])

            # Flip cell array?
            if len(cell_array.shape) == 2:
                flipped_array = cell_array[:, ::-1]
            else:
                flipped_array = cell_array[:, :, ::-1]

            # If necessary, average the data
            if averaged:
                final_array = np.zeros(netcdf_cell.shape[1:], dtype=np.float64)
                 # Compute the cell centered average
                standard_var_name = flopy_var_name.strip().replace('.', '_').replace(' ', '_').replace('-', '_').lower()
                self._average_flows(flipped_array, final_array, standard_var_name)
            else:
                final_array = np.copy(flipped_array)

            # Compare
            if not self._array_comp(final_array, netcdf_data):
                self._print_array_diff(final_array, netcdf_data, 'debug_cell_array_mf.txt', 'debug_cell_array_netcdf.txt')
                raise VerifyException('Data verification error: %s does not match at time %s.  Unmatched arrays dumped to debug_cell_array_*.txt.' % (cdf_var_name, str(time)))

    def _verify_head_data(self, objNCData, head_obj, fillvalue):
        time_values = head_obj.get_times()

        # Get head data from NetCDF file
        netcdf_head = objNCData.variables['head'][:]

        # Loop through all time values
        for i, time in enumerate(time_values):
            # Reconstruct head array from flopy object
            head_array = head_obj.get_data(totim=time)
            for f in self.fills:
                head_array[head_array == f] = fillvalue
            head_array[head_array <= -1e15] = fillvalue
            # Flip head array?
            head_array[:, :, :] = head_array[:, :, ::-1]

            # Compare to NetCDF
            if self._array_comp(head_array, netcdf_head[i, :, :, :].data) == False:
                self._print_array_diff(head_array, netcdf_head[i, :, :, :].data, 'debug_head_array_mf.txt', 'debug_head_array_netcdf')
                raise VerifyException('Data verification error: Heads does not match at time %d.  Unmatched arrays dumped to debug_cell_array_*.txt.' % time)

# General support methods #
    def _get_fill_values(self):

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

    def _get_coordinates(self):
        y, x, z = self.dis.get_node_coordinates()

        # Compute grid from top left corner, so switch Y back around
        y = y[::-1]

        # Compute grid from top left corner, so switch X back around
        x = x[::-1]

        # Convert Z to cartesian
        z = z[:, :, ::-1]

        # Convert distances to grid cell centers (from origin) to meters
        if self.grid_units == 'feet':
            x *= 0.3048
            y *= 0.3048
            z *= 0.3048

        logger.debug("Input origin point: {!s}, {!s}".format(self.grid_x, self.grid_y))

        # Support doing model grid location calculations in specified projected coordinate systems
        c_supported_crs_codes = [(0, 3999), (4399, 4462), (4484, 4489), (4491, 4554), (4568, 4589), (4652, 4656), (4766, 4800), (4855, 4880), (5014, 5016), (5105, 5130), (5167, 5188), (5253, 5259), (5292, 5316), (5343, 5349), (20004,32760)]
        supported_code = False
        if self.proj_epsg_code != self.c_epsg_code:
            if self.epsg_code == True:
                for crs_range in c_supported_crs_codes:
                    if crs_range[0] >= self.proj_epsg_code <= crs_range[1]:
                        supported_code = True
                        break
            else:
                supported_code = True

        # If in supported projected coordinate system: Calculate model grid points from northeast corner of
        # model grid using original coordinate system and then project into c_epsg_code
        if supported_code:
            print 'Calculating model grid points in original coordinate system.'
            # For supported codes X values increase to the east and Y values increase to the north
            # Calculate cell centers in projected coordinate system first and then convert to lat/long
            rotated_xsa = np.ndarray(0, dtype='float64')
            rotated_ysa = np.ndarray(0, dtype='float64')
            notrotated_xsa = np.ndarray(0, dtype='float64')
            notrotated_ysa = np.ndarray(0, dtype='float64')

            for x_val in x:
                # Build location matrices
                notrotated_xsa = np.append(notrotated_xsa, map(lambda val: self.grid_x + x_val, y))
                notrotated_ysa = np.append(notrotated_ysa, map(lambda val: self.grid_y - val, y))
                rotated_xsa = np.append(rotated_xsa, map(lambda val: self.grid_x + (x_val * math.cos(math.radians(self.grid_rotation)) - val * math.sin(math.radians(self.grid_rotation))), y))
                rotated_ysa = np.append(rotated_ysa, map(lambda val: self.grid_y - (x_val * math.sin(math.radians(self.grid_rotation)) + val * math.cos(math.radians(self.grid_rotation))), y))

            # Convert to lat/long
            notrotated_xsa, notrotated_ysa = self._transform_CRS_Matrix(notrotated_xsa, notrotated_ysa)
            rotated_xsa, rotated_ysa = self._transform_CRS_Matrix(rotated_xsa, rotated_ysa)

            # Convert to 2-D matrix
            xs_matrix = np.reshape(rotated_xsa, (len(x), len(y)))
            ys_matrix = np.reshape(rotated_ysa, (len(x), len(y)))

            # Assign to NetCDF variables
            self.netcdfdata["longitude"]  = rotated_xsa.reshape(self.dis.ncol, self.dis.nrow).T
            self.netcdfdata["latitude"] = rotated_ysa.reshape(self.dis.ncol, self.dis.nrow).T

        # Convert northeast corner of model grid to lat/long
        self.grid_x, self.grid_y = self._transform_CRS(self.grid_x, self.grid_y)

        # If not supported projected coordinate system: Project northeast corner of model grid into c_epsg_code
        # and calculate model grid point locations in c_epsg_code
        if not supported_code:
            print 'Calculating model grid points in output coordinate system.'

            # Convert to lat/long and then calculate cell centers directly in lat/long
            notrotated_xs = np.ndarray(0, dtype='float64')
            notrotated_ys = np.ndarray(0, dtype='float64')
            with LoggingTimer("Computing unrotated output grid", logger.info):
                upper = great_circle(distance=x, latitude=self.grid_y, longitude=self.grid_x, azimuth=90)
                for top_x, top_y in zip(upper["longitude"], upper["latitude"]):
                    # Compute the column points for each point in the upper row.
                    # Because this is a regular grid (rectangles), we can just rotate by 180 degrees plus the rotation angle.
                    row = great_circle(distance=y, latitude=top_y, longitude=top_x, azimuth=180)
                    notrotated_xs = np.append(notrotated_xs, row["longitude"])
                    notrotated_ys = np.append(notrotated_ys, row["latitude"])

            if self.grid_rotation != 0:
                with LoggingTimer("Computing rotated output grid", logger.info):
                    rotated_xs = np.ndarray(0, dtype='float64')
                    rotated_ys = np.ndarray(0, dtype='float64')
                    upper_rotated = great_circle(distance=x, latitude=self.grid_y, longitude=self.grid_x, azimuth=90+self.grid_rotation)
                    for top_x, top_y in zip(upper_rotated["longitude"], upper_rotated["latitude"]):
                        # Compute the column points for each point in the upper row.
                        # Because this is a regular grid (rectangles), we can just rotate by 180 degrees plus the rotation angle.
                        row = great_circle(distance=y, latitude=top_y, longitude=top_x, azimuth=180+self.grid_rotation)
                        rotated_xs = np.append(rotated_xs, row["longitude"])
                        rotated_ys = np.append(rotated_ys, row["latitude"])
            else:
                rotated_ys = notrotated_ys
                rotated_xs = notrotated_xs

            self.origin_x = self.grid_x
            self.origin_y = self.grid_y
            self.no_rotation_xs = notrotated_xs.reshape(self.dis.ncol, self.dis.nrow).T
            self.no_rotation_ys = notrotated_ys.reshape(self.dis.ncol, self.dis.nrow).T

            # Convert to 2-D matrix
            xs_matrix = np.reshape(rotated_xs, (len(x), len(y)))
            ys_matrix = np.reshape(rotated_ys, (len(x), len(y)))

            # Assign to NetCDF variables
            self.netcdfdata["longitude"]  = rotated_xs.reshape(self.dis.ncol, self.dis.nrow).T
            self.netcdfdata["latitude"] = rotated_ys.reshape(self.dis.ncol, self.dis.nrow).T

        self.netcdfdata["elevation"]  = z

    def compact_ibr(self, multi_array):
        compact_array = np.zeros(multi_array.shape)
        for layer in multi_array:
            row_idx = 0
            for row in layer:
                col_idx = 0
                for ibv in row:
                    if ibv != 0:
                        compact_array[0, row_idx, col_idx] = 1
                    col_idx += 1
                row_idx += 1
        return compact_array[0,:,::-1]

    def _transform_CRS_Matrix(self, x_matrix, y_matrix):
        x_matrix_tf = np.ndarray(0, dtype='float64')
        y_matrix_tf = np.ndarray(0, dtype='float64')
        for xval, yval in zip(x_matrix, y_matrix):
            tf = self._transform_CRS(xval, yval)
            x_matrix_tf = np.append(x_matrix_tf, tf[0])
            y_matrix_tf = np.append(y_matrix_tf, tf[1])
        return (x_matrix_tf, y_matrix_tf)

    def _transform_CRS(self, x_coord, y_coord):
        # CRS
        if (self.epsg_code == True and self.proj_epsg_code != self.c_epsg_code) or (self.epsg_code == False):
            # Transform origin_x and origin_y to correct projection
            if self.epsg_code:
                try:
                    grid_proj = Proj(init='epsg:{0!s}'.format(self.proj_epsg_code))
                except RuntimeError:
                    raise ValueError("Could not understand EPSG code '{!s}' from config file.".format(self.proj_epsg_code))
            else:
                try:
                    grid_proj = Proj(self.project_proj4str)
                except RuntimeError:
                    raise ValueError("Could not understand proj4 string'{!s}' from config file.".format(self.project_proj4str))

            wgs84_proj = Proj(init='epsg:%s' % self.c_epsg_code)
            transform_x, transform_y = transform(grid_proj, wgs84_proj, x_coord, y_coord)

            return transform_x, transform_y
        else:
            return x_coord, y_coord

    def _parse_config_file(self, config_file):
        if config_file is None:
            raise ValueError("No config file provided")

        config = ConfigParser.SafeConfigParser()
        try:
            config.read(config_file)
        except ConfigParser.ParsingError as e:
            raise ValueError("Bad configuration file.  Please check the contents. {!s}.".format(e.message))

        # Coordinate system
        self.epsg_code = False
        self.proj_epsg_code = None
        try:
            self.proj_epsg_code = config.get('space', 'epsg')
            self.epsg_code = True
        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as e:
            print 'No epsg code supplied.'
        if not self.epsg_code:
            try:
                self.project_proj4str = config.get('space', 'proj4str')
            except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as e:
                raise ValueError('No coorindate system supplied in configuration file.  Coordinate system must'
                                 'be supplied as an epsg code or a proj4 string.')

        try:
            self.grid_x        = config.getfloat('space', 'origin_x')
            self.grid_y        = config.getfloat('space', 'origin_y')
            self.grid_rotation = config.getfloat('space', 'rotation')

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
            self.base_date = config.get('time', 'base')

        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as e:
            raise ValueError('Configuration file is missing a required section: {0}'.format(e.message))

        # Output files (optional)
        try:
            self.cbud_file = os.path.join(os.path.dirname(config_file), config.get('output', 'cbud'))
        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as e:
            self.cbud_file = None

        try:
            self.head_file = os.path.join(os.path.dirname(config_file), config.get('output', 'head'))
        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as e:
            self.head_file = None

        try:
            self.head_type = config.get('output', 'headtype')
        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as e:
            self.head_type = 'binary'


        # Global attributes
        self.global_attributes = dict()
        if config.has_section('metadata'):
            self.global_attributes = { k : v for k, v in config.items('metadata') }

        return True
