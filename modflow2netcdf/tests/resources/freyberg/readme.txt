usage: mod2net.py [-h] -n NAME_FILE -c CONFIG_FILE [-o [OUTPUT_FILE]] [-v]
                  {plot,netcdf}

positional arguments:
  {plot,netcdf}

optional arguments:
  -h, --help            show this help message and exit
  -n NAME_FILE, --name_file NAME_FILE
                        Modflow Namefile (.nam) to load.
  -c CONFIG_FILE, --config_file CONFIG_FILE
                        Modflow2netcdf configiration file to load.
  -o [OUTPUT_FILE], --output_file [OUTPUT_FILE]
                        Output file (only used with the 'netcdf' action),
                        defaults to 'output.nc'
  -v, --verbose         Set output to be verbose, defaults to False




python ../../mod2net.py --name_file freyberg.nam -c freyberg.geo plot

python ../../mod2net.py --name_file freyberg.nam --config_file freyberg.geo netcdf

python ../../mod2net.py -n freyberg.nam -c freyberg.geo netcdf

python mod2net.py -n ./freyberg/modflowmodel/freyberg.nam -c ./freyberg/modflowmodel/freyberg.geo netcdf