# Modflow2NetCDF

### Version 1.0.0

## Introduction

Modflow2NetCDF is a tool for exporting geographic model data and results from MODFLOW 
projects to NetCDF format.  Modflow2NetCDF exports a MODFLOW project's model grid along with model 
cell elevation and location.  The geographic location and elevation of each model cell is calculated 
based on the model grid and additional information supplied in a configuration file. Model output 
data is exported for each cell from the head and cell by cell flow files.

## Installation

Modflow2NetCDF requires Python 2.7 or higher.  In addition, the following python libraries need to be installed prior to installation of Modflow2NetCDF.  Operating
system specific installation instructions are available below.

  * pyproj
  * pygc
  * NumPy
  * SciPy
  * matplotlib
  * HDF5
  * netCDF4
  * flopy

### Windows

1.  Install dependencies from http://www.lfd.uci.edu/~gohlke/pythonlibs/.  Take
special care to download the versions that match your version of Python.
  * pyproj
  * pygc
  * NumPy
  * SciPy
  * matplotlib
  * netCDF4

2.  Install 'flopy' (version 3) from https://github.com/modflowpy/flopy

### Linux (using pip)

1.  Install the following packages through your package management system or from source:
  * [HDF5 1.8.x](http://www.hdfgroup.org/HDF5/release/obtain5.html)
  * [netCDF 4.x](http://www.unidata.ucar.edu/downloads/netcdf/index.jsp) (with netCDF4/HDF5 support)

2.  Install the following libraries using 'pip install [library]'
  * pyproj (optional - only if you need to do a coordinate transformation)
  * pygc
  * numpy
  * scipy
  * matplotlib
  * netCDF4
  * flopy


### Linux (using conda)

1.  Install the following packages through your package management system or from source:
  * [HDF5 1.8.x](http://www.hdfgroup.org/HDF5/release/obtain5.html)
  * [netCDF 4.x](http://www.unidata.ucar.edu/downloads/netcdf/index.jsp) (with netCDF4/HDF5 support)
  * [PROJ.4](http://trac.osgeo.org/proj/)

2.  Install the following libraries using 'conda install [library]'
  * pyproj  (optional - only if you need to do a coordinate transformation)
  * pygc
  * numpy
  * scipy
  * matplotlib
  * netCDF4
  * flopy

## Running Modflow2NetCDF

The steps to run Modflow2NetCDF are:

1) Edit your output control file to generate the cell by cell flow and head output you wish to export
2) Run your model
3) Build a Modflow2NetCDF configuration file for you model (see example MODFLOW2NetCDF configuration file below)
4) Run mfnetcdf_cmdline.py (see example command line below)

For more information on steps 3 and 4 see the documentation and examples below.  In addition a tutorial ipython notebook is 
available with this distribution:

	docs\notebooks\MODFLOW NetCDF Visionalization.ipynb

## Documentation

MODFLOW2NetCDF supports a command line interface and a python library interface.  The command line interface
can be accessed by running the python script mfnetcdf_cmdline.py with the appropriate command line.
The library interface can be used by importing ModflowToNetCDF from modflow2netcdf.mfnetcdf.  Documentation
and examples for both methods are given below.

### Compatibility

MODFLOW2NetCDF's output NetCDF file can be displayed using the GODIVA2 Data Visualization option on a THREDDS data 
server.  GODIVA2 data visualization will only work properly when the data file has a consistent single time series.
Therefore head and cell by cell flow data must be saved during the same time intervals for these data to be displayed
correctly.  This can be accomplished by editing the MODFLOW output control file so that every where "SAVE HEAD" appears
in the file "SAVE BUDGET" also appears, and vica versa.

### Contents of NetCDF File

MODFLOW2NetCDF exports selected MODFLOW input and output data to a NetCDF file.  The contents of the exported NetCDF file
are as follows.

	* The length and width of the model cells in the MODFLOW grid (delc and delr arrays)
	* The size of the MODFLOW grid being exported (number of layers, rows, columns)
	* Three NetCDF dimension variables, describing the location of each model cell with x, y, and z coordinates
	* One NetCDF time dimension variable, containing all the times that data is present 

	* The following data sets may also be included, depending if a head and/or budget file are specified and the data available
	  in these files
		* Head values stored in the head file
		* Cell by cell flow values stored in the cell by cell flow file

### MODFLOW2NetCDF Configuration File

A MODFLOW2NetCDF configuration file needs to be built for each MODFLOW project exported into NetCDF format.  The 
configuration file contains information specific to the MODFLOW project being exported, including spatial and temporal 
information.  The configuration file is read using python's ConfigParser library, and consists of sections followed by 
"name: value" entries.

The configuration file is specified from the commandline (mfnetcdf_cmdline.py) with the -c [CONFIG_FILE] parameter.  The 
configuration file is specified through the library interface during initialization of the ModflowToNetCDF object with 
the config_file parameter.

Spatial data includes the grid projection used by the model, the location of the upper left grid point in 
the model's projection, the rotation of the grid from true north-south, and the units of measurement.  MODFLOW2NetCDF 
converts spatial data from your project's coordinate system into latitude and longitude 
coordinates on the WGS84 reference ellipsoid (EPSG 4326).  This allows all projects to be displayed
in a standard coordinate system.

Two example configuration files are shown below.  Configuration files are also available in two test projects
included with MODFLOW2NetCDF:

	tests\resources\freyberg\freyberg.geo

	tests\resources\carolina\carolina.geo

#### Configuration settings

The MODFLOW2NetCDF configuration file has five sections, general, space, time, output, and metadata.  The
sections and descriptions of required and optional entries in each section are documented below.

##### General section
###### precision : single or double
	Precision that the MODFLOW output file was saved in
	
##### Space section
###### crs : [integer or string]
	EPSG code or pyproj4 string of the grid projection used by the model
###### origin_x : [float]
	X coordinate location of the upper left corner of the model grid in the projection used by the model
###### origin_y : [float]
	Y coordinate location of the upper left corner of the model grid in the projection used by the model
###### rotation : [float]
	Clockwise rotation in degrees of the model grid from true north
###### units : [meters or feet]
	Model's units of measurement
	
##### Time section
###### units: [string]
	Time units specified in the NetCDF file.  
		Example: 'days', 'hours', or 'minutes'
###### base: [string]
	Base date of when the model started.  UTC is assumed if no timezone information is specified.
		Example: '2006-06-01 00:00:00'

##### Output section

###### head: [file path] (optional)
	Path to the Head output file relative to configuration file.
		Example: 'output\mymodelrun.hds'

###### cbud: [file path] (optional)
	Path to the CellBudget output file relative to configuration file.
		Example: 'output\mymodelrun.cbb'

###### headtype: [binary or formatted] (optional)
	Type of head file, binary or formatted.
		Example: 'binary'

##### Metadata section

###### [key]: [value]
	Each key/value in the 'metadata' block will be added as a global attribute in the NetCDF4 file
		Example: 'creator:  modflow2netcdf'

### Command line interface

The 'mfnetcdf_cmdline.py' command line python script provides a simple interface to the MODFLOW2NetCDF library
that allows for the creation of a NetCDF file.  Usage of this interface requires only a basic understanding of 
using command line interfaces and does not require any python programming language knowledge.  The 'mfnetcdf_cmdline.py'
script takes the following command line parameters:

	-n NAME_FILE				
		The path to a MODFLOW namefile
	-c CONFIG_FILE			
		The path to a MODFLOW2NetCDF configuration file
	[-o OUTPUT_FILE]		
		The path to the MODFLOW2NetCDF output file, default output file is output.nc 
	[-v]								
		Run in verbose mode, defaults to False
	[-vf]								
		Verify NetCDF file, defaults to False
	
### Using Modflow2NetCDF as a library 

The Modflow2NetCDF library was designed to export data and results from a MODFLOW project to a NetCDF formatted file. The 
Modflow2NetCDF library contains a single class, ModflowToNetCDF.

#### ModflowToNetCDF Class

The ModflowToNetCDF class supports exporting MODFLOW data to a NetCDF file and plotting MODFLOW data.  This class reads a MODFLOW 
project's input and output files using flopy.  It also reads a Modflow2NetCDF configuration file that must be created for each 
MODFLOW project being exported into NetCDF format.

##### Parameters
	
		namfilename : string
			MODFLOW project name file including full or relative path
		config_file : string
			ModflowToNetCDF config file including full or relative path
		version : string (optional)
			Version of MODFLOW being used 
		  	default value = 'mf2k'
		exe_name : string (optional)
			Modflow executable name 
				default value = 'mf2005.exe'
		verbose : boolean (optional)
			Run ModflowtoNetCDF in verbose mode 
				default value = False
		model_ws : string (optional)
			Path to model input files including full or relative path 
				defualt_value = None (current path)
		
##### Methods
		
###### to_plot :  Creates a plot.
       
        Parameters
        ----------
				variable : string (optional)
					Variable name of data to plot 
						default value = None (plot model surface)
				level : integer (optional)
					Value representing the layer to plot 
						default value = 0 (top layer)
				time : integer (optional)
					The time the data will be plotted at 
						default value = 0
				colormap : Colormap (optional)
					Colormap to use for the plot 
						default value = matplotlib.cm.Reds
					
###### save_netcdf : Creates a netcdf file.

        Parameters
        ----------
        output_file : string
        	File name including full file path of output netcdf file
        verify : boolean (optional)
        	Netcdf file verification mode.  When set to True data in the netcdf output file is verified.
      			default value = False     
      			 
        Returns
        ----------
        netcdf4.dataset
            dataset containing the contents of the outputted netcdf file
 

## Examples

### Modflow2NetCDF Command Line Example

Modflow2NetCDF can be run from the command line with three required command line switches 
identifing the name file, Modflow2NetCDF config file, and output file:

	mfnetcdf_cmdline.py -n "input\modflow_proj.nam" -c "input\modflow_proj.geo" -o "output\netcdf_file.nc"

See the command line example in docs\examples\commandline.
	
### Modflow2NetCDF Configuration File Examples

The Modflow2NetCDF configuration file contains project specific information about a specific MODFLOW project which is
used to properly export MODFLOW data into NetCDF format.  The following are examples of Modflow2NetCDF configuration
files.  These examples contain project specific information that most likely will need to be changed for each MODFLOW 
project. 

#### Example WGS84 Configuration File
```ini
[general]
precision:  double    ; If the model was run with single or double precision. 'single' or 'double'.

[space]
crs:      4326        ; EPSG code of grid projection used by the model
origin_x: -72         ; Longitude of upper left grid point (origin)
origin_y: 41.44       ; Latitude of upper left grid point (origin)
rotation: 45          ; True north based grid rotation angle (clockwise from true north)
units:    m           ; Units of measurment in output ('meters', 'm', 'feet', 'ft', or 'f')

[time]
units:    days        ; Time units in output
base:     2006-06-01 00:00:00    ; Assumed UTC if no timezone information is specified

[output]
head:     mymodelrun.hds  ; Optional. Path to the Head output file (relative to config file).
cbud:     mymodelrun.cbb  ; Optional. Path to the CellBudget output file (relative to config file).

[metadata]            ; Each key/value in the 'metadata' block will be added as a global attribute in the NetCDF4 file
id:       my_model_id
creator:  modflow2netcdf

```

#### Example Web-mercator Configuration File
```ini
[general]
precision:  single    ; If the model was run with single or double precision. 'single' or 'double'.

[space]
crs:      3857        ; EPSG code of grid projection used by the model
origin_x: -8012405.88 ; Longitude of upper grid left point (origin)
origin_y: 5078408.56  ; Latitude of uppper grid left point (origin)
rotation: 0           ; True north based grid rotation angle (clockwise from true north)
units:    ft          ; Units of measurment in output ('meters', 'm', 'feet', 'ft', or 'f')

[time]
units:    days        ; Time units in output
base:     1992-01-06 06:00:00 -0500

[output]  ; Not needed, the default extensions for the Head (.bhd) and CellBudget (.bud) output files will be assumed.

[metadata]            ; Each key/value in the 'metadata' block will be added as a global attribute in the NetCDF4 file
id:       my_model_id
creator:  modflow2netcdf
```

#### Example Python Code Accessing the MODFLOW2NetCDF Library

MODFLOW2NetCDF provides a python library interface.  An example python script using this interface can be found in:

	docs\examples\script\mfnetcdf_example.py
	
This example script exports MODFLOW data from the Freyberg test model to a NetCDF4 file named freyberg.nc.  The script 
makes use of freyberg.geo, a Modflow2NetCDF configuration file written for the Freyberg project.

## Testing

Unit tests can be found in the 'tests' folder.  Unit tests are executed by running the unit test code in 'mfnetcdftest.py'.
Unit tests are implemented using python's unittest library.  Detailed results of the unit tests are written to 
TestProjectLog.txt.

Unit tests run two types of tests, internal data verification tests and external data verification tests.   

Internal data verification tests create a NetCDF file and then read the NetCDF file using the NetCDF4 interface.  The
contents of the NetCDF file is then verified against the MODFLOW project's data, accessed using FLOPY.  Internal data verification
tests are designed to find errors during the conversion and export of the MODFLOW data to NetCDF format.

External data verification tests verify the contents of a NetCDF file against externally produced (produced without the use
of Modflow2NetCDF or any of its dependent libraries) expected results.  While external data verification tests are more rigorous,
creating these tests is also more time intensive.  Therefore there are only a limited number of external data verification tests
available for Modflow2NetCDF.