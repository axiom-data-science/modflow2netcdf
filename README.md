# Netcdf_cmdline

### Version 1.0.0

## Introduction

Netcdf_cmdline is a command line interface to flopy's netcdf export feature.  Netcdf_cmdline is used
to export geographic model data and results from MODFLOW projects to NetCDF format.  Netcdf_cmdline 
uses flopy to export a MODFLOW project's model grid along with model cell elevation and location.  The 
geographic location and elevation of each model cell is calculated based on the model grid and additional 
information supplied in a configuration file.  Model input and output data are exported for each cell from 
the modflow input and output files.

Netcdf_cmdline uses the Flopy libraries to access MODFLOW project data, and supports the versions
of MODFLOW supported by Flopy, including MODFLOW-2000, MODFLOW-2005, MODFLOW-NWT, and MODFLOW-USG.

## Installation

Netcdf_cmdline requires Python 2.7 or higher.  In addition, the following python libraries need to be installed prior to using netcdf_cmdline.  Operating
system specific installation instructions are available below.
  * pyproj
  * NumPy
  * HDF5
  * netCDF4
  * flopy

### Windows

1.  Install dependencies from http://www.lfd.uci.edu/~gohlke/pythonlibs/.  Take
special care to download the versions that match your version of Python.
  * pyproj
  * NumPy
  * HDF5
  * netCDF4
  * flopy

2.  Install 'flopy' (version 3) from https://github.com/modflowpy/flopy

### Linux (using pip)

1.  Install the following packages through your package management system or from source:
  * [HDF5 1.8.x](http://www.hdfgroup.org/HDF5/release/obtain5.html)
  * [netCDF 4.x](http://www.unidata.ucar.edu/downloads/netcdf/index.jsp) (with netCDF4/HDF5 support)

2.  Install the following libraries using 'pip install [library]'
  * pyproj
  * NumPy
  * HDF5
  * netCDF4
  * flopy


### Linux (using conda)

1.  Install the following packages through your package management system or from source:
  * [HDF5 1.8.x](http://www.hdfgroup.org/HDF5/release/obtain5.html)
  * [netCDF 4.x](http://www.unidata.ucar.edu/downloads/netcdf/index.jsp) (with netCDF4/HDF5 support)
  * [PROJ.4](http://trac.osgeo.org/proj/)

2.  Install the following libraries using 'conda install [library]'
  * pyproj
  * NumPy
  * HDF5
  * netCDF4
  * flopy

## Running Netcdf_cmdline

The steps to run Netcdf_cmdline are:

1) Add geographic location information to your project's MODFLOW name file

2) Edit your project's MODFLOW output control file to generate the cell budget and head output you wish to export

3) Run your MODFLOW model

4) Run Netcdf_cmdline.py (see example command line in the "Examples" section)

## Documentation

Netcdf_cmdline is a command line interface and a netcdf export feature included in flopy's python library interface.  
The command line interface is accessed by running the python script Netcdf_cmdline.py with the appropriate command line.

### Compatibility

Netcdf_cmdline generates twp output NetCDF file that can be displayed using the GODIVA2 Data Visualization option on a 
THREDDS data server.  GODIVA2 data visualization will only work properly for data files that have a consistent single 
time series.  Therefore head and cell budget data must be saved during the same time intervals for these data to be displayed
correctly.  This can be accomplished by editing the MODFLOW output control file so that for every stress period and time
step that "SAVE HEAD" appears, "SAVE BUDGET" also appears, and vica versa.

### Command line interface

The 'Netcdf_cmdline.py' python script provides a simple command line interface to flopy's netcdf export library.  This 
interface can be used to export data from a MODFLOW project to a NetCDF file.  Usage of this interface requires only a 
basic understanding of command line interfaces and does not require any python programming language knowledge.  The 
'Netcdf_cmdline.py' script takes the following command line parameters:

	-n NAME_FILE				
		The path to a MODFLOW namefile
	-ni FILE_PATH
		The path to the NetCDF output file that will store model input data
	-no FILE_PATH
		The path to the NetCDF output file that will store model output data
	-p PRECISION ('single' or 'double')
		Precision of modflow output files, single or double.  

### Geographic location information

Geographic location information must be added to the beginning of your project's MODFLOW name file.  Geographic location 
information includes x and y coordinates of the northwest (upper left) corner of the model grid (xul and yul), the model 
grid's rotation, and the proj4 coordinate system of your model grid (proj4_str).  Geographic location information is added 
to a single commented line in the MODFLOW name file using the following format:

\#xul:[X Cooordinate], yul:[Y Coordinate], rotation:[Rotation Angle], proj4_str:[Proj4 String]
		
## Examples

### Netcdf_cmdline Command Line Example

Netcdf_cmdline can be run using four required command line switches identifing the name file, the NetCDF output file for
MODFLOW input data, the NetCDF output file for MODFLOW output data, and the precision of the MODFLOW output files:

	mfnetcdf_cmdline.py -n "input\modflow_proj.nam" -ni "ouput\netcdf_file_in.nc" -no "output\netcdf_file_out.nc" -p double

### Geographic information

Below is an example of geographic coorindate information added to the top of a MODFLOW name file.  The example coordinate 
information is for a model grid with x and y coordinate locations (northwest corner) of -83.747389 and 32.917647, 
grid rotation of -42.95 degrees, using the latitute/longitude WGS84 coorindates.

	\# MF2K NAME file
	\#xul:-83.747389, yul:32.917647, rotation:-42.95, proj4_str:+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs

## NetCDF Output File Contents

The NetCDF output file contains location based data saved with latitute, longitude, and elevation coordinates based on EPSG 
code 4326.  Data is layed out on a grid with spatial dimension variables (x, y, and layer).  In addition, a time dimension (t) 
is used to specify the times used in the MODFLOW head and cell budget output files.  The NetCDF output file contains the 
following variables.

##### Dimensions
	
	x
	y
	layer
	time
		
##### MODFLOW Input Variables Stored

MODFLOW input variables stored in the NetCDF output file include the following.  In addition, many package specific
input variables are also stored.
	
	crs 
		Coordinate system used.  EPSG code: 4326
	latitude
		2-D array of latitude values for each model cell in a single model layer
	longitude
		2-D array of longitude values for each model cell in a single model layer
	time
		1-D array of times for output data (head values and cell budget values)
	elevation
		3-D array of elevation values for each model cell
	layer
		1-D array of model layers
	delc
		1-D array of cell widths along model columns
	delr
		1-D array of cell widths along model rows
	layer thickness
		3-D array of layer thickness
	starting heads
		3-D array of starting heads
	ibound
		3-D ibound array
		
##### MODFLOW Output Variables Stored

MODFLOW output variables stored in the NetCDF output file are read from the output head file and cell budget
file.  They include the following:
	
Head file
	head
		4-D array of head values (time, layer, x, y)

Cell by cell Budget File
	constant_head
	flow_right_face
	flow_front_face 
	flow_lower_face 
	wells
	drains
	river_leakage
	head_dep_bounds
	recharge
	specified_flows
	stream_leakage
	et_segments
	mnw
	storage