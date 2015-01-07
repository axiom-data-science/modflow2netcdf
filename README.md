## Modflow2NetCDF (mod2net)

### Installation

##### Windows

1.  Install dependencies from http://www.lfd.uci.edu/~gohlke/pythonlibs/.  Take
special care to download the versions that match your version of Python.
  * pyproj
  * NumPy
  * SciPy
  * matplotlib
  * netCDF4

2.  Install `flopy` (version 3) from https://github.com/modflowpy/flopy

3.  Install the following libraries using `pip install [library]`
  * pytz
  * pygc
  * python-dateutil

##### Linux (using pip)

1.  Install the following packages through your package management system or from source:
  * [HDF5 1.8.x](http://www.hdfgroup.org/HDF5/release/obtain5.html)
  * [netCDF 4.x](http://www.unidata.ucar.edu/downloads/netcdf/index.jsp) (with netCDF4/HDF5 support)
  * [PROJ.4](http://trac.osgeo.org/proj/)

2.  Install the following libraries using `pip install [library]`
  * pyproj
  * numpy
  * scipy
  * matplotlib
  * netCDF4
  * flopy
  * pygc
  * python-dateutil
  * pytz


##### Linux (using conda)

1.  Install the following packages through your package management system or from source:
  * [HDF5 1.8.x](http://www.hdfgroup.org/HDF5/release/obtain5.html)
  * [netCDF 4.x](http://www.unidata.ucar.edu/downloads/netcdf/index.jsp) (with netCDF4/HDF5 support)
  * [PROJ.4](http://trac.osgeo.org/proj/)

2.  Install the following libraries using `conda install [library]`
  * pyproj
  * numpy
  * scipy
  * matplotlib
  * netCDF4
  * flopy
  * pygc
  * python-dateutil
  * pytz


### Running

#### Inputs

The `mod2net` command takes in two parameters, the path to a Modflow namefile
(.nam) and the path to a `mod2net` [configuration file](#configuration-file) file.

###### Examples


#### Configuration File

A `mod2net` configuration file should be formated as so:


###### WGS84 Configuration
```ini
[general]
precision:  double    ; If the model was run with single or double precision. 'single' or 'double'.

[space]
crs:      4326        ; EPSG code of grid projection
origin_x: -72         ; Longitude of upper grid left point (origin)
origin_y: 41.44       ; Latitude of uppper grid left point (origin)
rotation: 45          ; True north based grid rotation angle (clockwise from true north)
units:    m           ; Units of measurment in output ('meters', 'm', 'feet', 'ft', or 'f')

[time]
units:    days
base:     2006-06-01 00:00:00    ; Assumed UTC if no timezone information is specified

[output]
head:     mymodelrun.hds  ; Optional. Path to the Head output file (relative to config file).
cbud:     mymodelrun.cbb  ; Optional. Path to the CellBudget output file (relative to config file).

```

###### Web-mercator Configuration
```ini
[general]
precision:  single    ; If the model was run with single or double precision. 'single' or 'double'.

[space]
crs:      3857        ; EPSG code of grid projection
origin_x: -8012405.88 ; Longitude of upper grid left point (origin)
origin_y: 5078408.56  ; Latitude of uppper grid left point (origin)
rotation: 0           ; True north based grid rotation angle (clockwise from true north)
units:    ft          ; Units of measurment in output ('meters', 'm', 'feet', 'ft', or 'f')

[time]
units:    days
base:     1992-01-06 06:00:00 -0500

[output]  ; Not needed, the default extensions for the Head (.bhd) and CellBudget (.bud) output files will be assumed.
```


#### Testing

Testing is done through `pytest` excecuted in the root of the source code.

```bash
python -m pytest
```

To run the tests for the Denver Basin, the output files must be present in the `modflow2netcdf/tests/resources/colorado` directory.  Output can be downloaded from http://pubs.usgs.gov/pp/1770/.  Copy all files from `C-Transient-final/Calibrated_model/` into the `modflow2netcdf/tests/resources/colorado` folder.
