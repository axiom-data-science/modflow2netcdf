# Modflow2NetCDF

### Version 1.0.0

## Introduction

Modflow2NetCDF is a tool for exporting geographic model data and results from MODFLOW 
projects to NetCDF format.  Modflow2NetCDF exports a MODFLOW project's model grid along with model 
cell elevation and location.  The geographic location and elevation of each model cell is calculated 
based on the model grid and additional information supplied in a configuration file.  Model output 
data is exported for each cell from the head and cell budget files.

Modflow2NetCDF uses the Flopy libraries to access MODFLOW project data, and supports the versions
of MODFLOW supported by Flopy, including MODFLOW-2000, MODFLOW-2005, MODFLOW-NWT, and MODFLOW-USG.

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

1) Edit your MODFLOW output control file to generate the cell budget and head output you wish to export

2) Run your MODFLOW model

3) Build a Modflow2NetCDF configuration file for you model (see example MODFLOW2NetCDF configuration file in the "Examples" section)

4) Run mfnetcdf_cmdline.py (see example command line in the "Examples" section)

For more information on steps 3 and 4 see the documentation and examples below.  In addition a tutorial ipython notebook is 
available with this distribution:

	docs\notebooks\MODFLOW NetCDF Visionalization.ipynb

## Documentation

MODFLOW2NetCDF supports a command line interface and a python library interface.  The command line interface
is accessed by running the python script mfnetcdf_cmdline.py with the appropriate command line.
To use the library interface import ModflowToNetCDF from modflow2netcdf.mfnetcdf.  Documentation
and examples for both methods are given below.

### Compatibility

MODFLOW2NetCDF's output NetCDF file can be displayed using the GODIVA2 Data Visualization option on a THREDDS data 
server.  GODIVA2 data visualization will only work properly for data files that have a consistent single time series.
Therefore head and cell budget data must be saved during the same time intervals for these data to be displayed
correctly.  This can be accomplished by editing the MODFLOW output control file so that for every stress period and time
step that "SAVE HEAD" appears, "SAVE BUDGET" also appears, and vica versa.

### MODFLOW2NetCDF Configuration File

A MODFLOW2NetCDF configuration file needs to be built for each MODFLOW project exported into NetCDF format.  The 
configuration file contains information specific to the MODFLOW project being exported, including spatial and temporal 
information.  Python's ConfigParser library is used to read the configuration file, which consists of sections followed by 
"name: value" entries.

The configuration file is specified from the MODFLOW2NetCDF commandline (mfnetcdf_cmdline.py) with the -c [CONFIG_FILE] parameter.  
The MODFLOW2NetCDF library requires the location of the configuration file be specified during initialization of the 
ModflowToNetCDF class.

The configuration file contains spatial data which includes the grid projection used by the model, the location of the 
upper left grid point in the model's projection, the rotation of the grid from true north-south, and the units of measurement.  
MODFLOW2NetCDF converts spatial data from your project's coordinate system into latitude and longitude coordinates on the 
WGS84 reference ellipsoid (EPSG 4326).  This allows all projects to be displayed in a standard coordinate system.

Two example configuration files are shown in the "Examples" section.  Configuration files are also available with two test 
projects included with MODFLOW2NetCDF:

	tests\resources\freyberg\freyberg.geo

	tests\resources\carolina\carolina.geo

#### Configuration settings

The MODFLOW2NetCDF configuration file has five sections, general, space, time, output, and metadata.  The
sections and descriptions of required and optional entries in each section are documented below.

##### General section
###### precision : [single or double]
	Precision of the MODFLOW output file
	
##### Space section
###### crs : [integer or string]
	EPSG code or pyproj4 string of the grid projection used by the MODFLOW model
###### origin_x : [float]
	X coordinate location of the upper left corner of the model grid in the projection used by the MODFLOW model
###### origin_y : [float]
	Y coordinate location of the upper left corner of the model grid in the projection used by the MODFLOW model
###### rotation : [float]
	Clockwise rotation of the model grid in degrees from true north
###### units : [meters or feet]
	Model's units of measurement
	
##### Time section
###### units: [string]
	Time units specified in the NetCDF file.  
		Example: 'days', 'hours', or 'minutes'
###### base: [string]
	Base date when the model started.  UTC is assumed if no timezone information is specified.
		Example: '2006-06-01 00:00:00'

##### Output section

###### head (optional): [file_path] 
	Path to the head output file relative to configuration file.
		Example: 'output\mymodelrun.hds'

###### cbud (optional): [file_path]
	Path to the cellBudget output file relative to configuration file.
		Example: 'output\mymodelrun.cbb'

###### headtype (optional): [binary/formatted]
	Type of head file (binary or formatted).
		Example: 'binary'

##### Metadata section

###### [key]: [value]
	Each key/value in the 'metadata' block will be added as a global attribute in the NetCDF4 file
		Example: 'creator:  modflow2netcdf'

### Command line interface

The 'mfnetcdf_cmdline.py' python script provides a simple command line interface to the MODFLOW2NetCDF library.  This 
interface can be used to export data from a MODFLOW project to a NetCDF file.  Usage of this interface requires only a 
basic understanding of command line interfaces and does not require any python programming language knowledge.  The 
'mfnetcdf_cmdline.py' script takes the following command line parameters:

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

The ModflowToNetCDF class contains two methods, one method exports MODFLOW data to a NetCDF file and one method plots 
MODFLOW data.  This class reads data from a MODFLOW project's input and output files using flopy.  This class requires a 
Modflow2NetCDF configuration file that must be created for each MODFLOW project that is exported into NetCDF format.

##### Parameters
	
		namfilename : string
			MODFLOW project name file including full or relative path
		config_file : string
			ModflowToNetCDF config file including full or relative path
		version : string (optional)
			Version of MODFLOW project 
		  	default value = 'mf2k'
		exe_name : string (optional)
			Modflow executable name 
				default value = 'mf2005.exe'
		verbose : boolean (optional)
			Run ModflowtoNetCDF in verbose mode 
				default value = False
		model_ws : string (optional)
			Full or relative path to model input files 
				defualt_value = None (current path)
		
##### Methods
		
###### to_plot :  Creates a plot.
       
        Parameters
        ----------
		variable : string (optional)
			Data variable name to plot 
				default value = None (plot model surface)
		level : integer (optional)
			Number of model layer to plot 
				default value = 0 (top layer)
		time : integer (optional)
			Date/time to be plotted
				default value = 0
		colormap : Colormap (optional)
			Colormap to use 
				default value = matplotlib.cm.Reds
					
###### save_netcdf : Creates a netcdf file.

        Parameters
        ----------
        output_file : string
        	File name including full file path of output netcdf file
        verify : boolean (optional)
        	Netcdf file verification mode.  When set to True data in the netcdf output file are verified.
      			default value = False     
      			 
        Returns
        ----------
        netcdf4.dataset
            Dataset containing the contents of the outputted netcdf file
 

## Examples

### Modflow2NetCDF Command Line Example

Modflow2NetCDF can be run from the command line using three required command line switches 
identifing the name file, Modflow2NetCDF config file, and output file:

	mfnetcdf_cmdline.py -n "input\modflow_proj.nam" -c "input\modflow_proj.geo" -o "output\netcdf_file.nc"

See the command line example in docs\examples\commandline.
	
### Modflow2NetCDF Configuration File Examples

The Modflow2NetCDF configuration file contains project specific information about a specific MODFLOW project which is
used to properly export MODFLOW data into NetCDF format.  The following are examples of Modflow2NetCDF configuration
files.  These examples contain project specific information that most likely will need to be customized for each MODFLOW 
project. 

#### Example WGS84 Configuration File

[general]
precision:  double    ; The precision of the MODFLOW model, single or double precision. ('single' or 'double').

[space]
crs:      4326        ; EPSG code of grid projection used by the model
origin_x: -72         ; Longitude of upper left grid point (origin)
origin_y: 41.44       ; Latitude of upper left grid point (origin)
rotation: 45          ; True north based grid rotation angle (clockwise from true north)
units:    m           ; Units of measurment ('meters', 'm', 'feet', 'ft', or 'f')

[time]
units:    days        ; Time units in output
base:     2006-06-01 00:00:00    ; Assumed UTC if no timezone information is specified

[output]
head:     mymodelrun.hds  ; Optional. Path to the head output file (relative to config file).
cbud:     mymodelrun.cbb  ; Optional. Path to the cell budget output file (relative to config file).

[metadata]            ; Each key/value in the 'metadata' block will be added as a global attribute in the NetCDF4 file
id:       my_model_id
creator:  modflow2netcdf

```

#### Example Web-mercator Configuration File
```ini
[general]
precision:  single    ; The precision of the MODFLOW model, single or double precision. ('single' or 'double').

[space]
crs:      3857        ; EPSG code of grid projection used by the model
origin_x: -8012405.88 ; Longitude of upper grid left point (origin)
origin_y: 5078408.56  ; Latitude of uppper grid left point (origin)
rotation: 0           ; True north based grid rotation angle (clockwise from true north)
units:    ft          ; Units of measurment ('meters', 'm', 'feet', 'ft', or 'f')

[time]
units:    days        ; Time units in output
base:     1992-01-06 06:00:00 -0500

[output]  ; Not needed, the default extensions for the head (.bhd) and cell budget (.bud) output files will be assumed.

[metadata]            ; Each key/value in the 'metadata' block will be added as a global attribute in the NetCDF4 file
id:       my_model_id
creator:  modflow2netcdf
```

#### Example Python Code Accessing the MODFLOW2NetCDF Library

MODFLOW2NetCDF provides a python library interface.  An example python script using this interface can be found here:

	docs\examples\script\mfnetcdf_example.py
	
This example script exports MODFLOW data from the Freyberg test model to a NetCDF4 file named freyberg.nc.  The script 
makes use of freyberg.geo, a Modflow2NetCDF configuration file written for the Freyberg project.

## NetCDF Output File Contents

The NetCDF output file that Modflow2NetCDF creates contains location based data saved with latitute, longitude, and elevation
coordinates based on EPSG code 4326 (semi-major axis = 6378137.0, inverse flattening = 298.257223563).  Data is layed out on a 
grid with spatial dimension variables (x, y, and layer).  In addition, a time dimension (t) is used to specify the times used 
in the MODFLOW head and cell budget output files.  The NetCDF output file contains the following variables.

##### Dimensions
	
	x
	y
	layer
	time
		
##### General Variables
	
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
	VerticalTransform
		
##### Head variables (when provided in head file)
	
	head
		4-D array of head values (time, layer, x, y)

##### Cell budget variables (when provided in cell budget file)
	
	Cell budget variables each are written as 4-D arrays (time, layer, x, y) to the 
	NetCDF output file.  These variables may include:
		
		constant_head
		flow_right_face_centered 
		flow_right_face
		flow_front_face_centered
		flow_front_face 
		flow_lower_face_centered
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
		
	All variables except the flow_xxxx_face_centered variables are extracted directly from the cell budget file.  The flow_xxxx_face_centered
	variables are calculated 4-D arrays based on the cooresponding flow_xxxx_face variables extracted from the cell budget file.  The centered
	values are calculated as the average flow across two adjacent (parallel) faces.  For example, flow_right_face_centered is calculated for cell
	(i,j,k) as the average of the flow across the right face of (i,j,k) (the face shared by (i,j,k) and (i,j+1,k)) and the flow across the right
	face of (i,j-1,k) (the face shared by (i,j-1,k) and (i,j,k). 

## Testing

Unit tests can be found in the 'tests' folder.  Unit tests are executed by running the unit test code in 'mfnetcdftest.py'.
Unit tests are implemented using python's unittest library.  Detailed results of the unit tests are written to 
TestProjectLog.txt.

Unit tests run two types of tests, internal data verification tests and external data verification tests.   

Internal data verification tests occur after creation of a NetCDF file.  The test read the NetCDF file using the NetCDF4 interface.  The
contents of the NetCDF file is then verified against the MODFLOW project's data, accessed using the FLOPY libraries.  Internal data verification
tests are designed to find errors during the conversion and export of the MODFLOW data to NetCDF format.

External data verification tests verify the contents of a NetCDF file against externally produced (produced without the use
of Modflow2NetCDF or any of its dependent libraries) expected results.  External data verification tests are more rigorous, but can only
be run on projects where expected results have been previously generated by an unrelated process.  Therefore there are only a limited number 
of external data verification tests available for Modflow2NetCDF.