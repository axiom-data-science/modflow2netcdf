import os
import logging
import unittest

from modflow2netcdf.output import ModflowOutput
import netCDF4

from modflow2netcdf import logger
logger.level = logging.DEBUG
logger.addHandler(logging.StreamHandler())


class TestOutput(unittest.TestCase):

    def setUp(self):
        self.resource      = os.path.join(os.path.dirname(__file__), "resources", "l1a2k.nam")
        self.wgs84_geofile = os.path.join(os.path.dirname(__file__), "resources", "wgs84.geo")
        self.web_geofile   = os.path.join(os.path.dirname(__file__), "resources", "web.geo")
        self.bad_geofile   = os.path.join(os.path.dirname(__file__), "resources", "bad.geo")
        self.output_file   = os.path.join(os.path.dirname(__file__), "resources", "output.nc")

    def test_no_geofile(self):
        with self.assertRaises(ValueError):
            ModflowOutput(self.resource, geofile=None, exe_name="mf2005", verbose=False)

    def test_bad_geofile(self):
        with self.assertRaises(ValueError):
            ModflowOutput(self.resource, geofile=self.bad_geofile, exe_name="mf2005", verbose=False)

    def test_load(self):
        mf = ModflowOutput(self.resource, geofile=self.wgs84_geofile, exe_name="mf2005", verbose=True)
        assert mf is not None

    def test_plot(self):
        mf = ModflowOutput(self.resource, geofile=self.wgs84_geofile, exe_name="mf2005", verbose=True)
        assert mf is not None
        mf.to_plot()

    def test_web_geofile(self):
        mf = ModflowOutput(self.resource, geofile=self.web_geofile, exe_name="mf2005", verbose=False)
        assert mf is not None
        mf.to_netcdf(output_file=self.output_file)
        nc = netCDF4.Dataset(self.output_file)
        assert nc is not None

    def test_convert(self):
        mf = ModflowOutput(self.resource, geofile=self.wgs84_geofile, exe_name="mf2005", verbose=False)
        assert mf is not None
        mf.to_netcdf(output_file=self.output_file)
        nc = netCDF4.Dataset(self.output_file)
        assert nc is not None

    def test_colorado_netcdf(self):
        colorado_nam = os.path.join(os.path.dirname(__file__), "resources", "colorado", "mod16_ssfix_wel2.nam")
        colorado_geo = os.path.join(os.path.dirname(__file__), "resources", "colorado", "colorado.geo")
        mf = ModflowOutput(colorado_nam, geofile=colorado_geo, exe_name="mf2005", verbose=True)
        assert mf is not None
        mf.to_netcdf(output_file=self.output_file)
        nc = netCDF4.Dataset(self.output_file)
        assert nc is not None
