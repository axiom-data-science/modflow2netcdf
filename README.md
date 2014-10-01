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

2.  Install `flopy` from http://code.google.com/p/flopy/
http://code.google.com/p/flopy/


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

### Running

#### Inputs

The `mod2net` command takes in two parameters, the path to a Modflow namefile
(.nam) and the path to a `mod2net` [configuration file](#configuration-file) file.

###### Examples


#### Configuration File

A `mod2net` configuration file should be formated as so:


###### WGS84 Configuration
```python
4326        # EPSG code of grid projection
-72         # Longitude of upper grid left point (origin)
41.44       # Latitude of uppper grid left point (origin)
45          # True north based grid rotation angle (clockwise from true north)
m           # Units of measurment in output ('meters', 'm', 'feet', or 'f')
```

###### Web-mercator Configuration
```python
3857        # EPSG code of grid projection
-8012405.88 # Longitude of upper grid left point (origin)
5078408.56  # Latitude of uppper grid left point (origin)
0           # True north based grid rotation angle (clockwise from true north)
ft          # Units of measurment in output ('meters', 'm', 'feet', 'ft', or 'f')
```

