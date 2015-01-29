import os
import logging
import unittest

from modflow2netcdf.output import ModflowOutput
import netCDF4

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ioff()

from modflow2netcdf import logger
logger.level = logging.DEBUG
logger.addHandler(logging.StreamHandler())


class TestOutput(unittest.TestCase):

    def setUp(self):
        self.resource          = os.path.join(os.path.dirname(__file__), "resources", "l1a2k.nam")
        self.wgs84_config_file = os.path.join(os.path.dirname(__file__), "resources", "wgs84.geo")
        self.web_config_file   = os.path.join(os.path.dirname(__file__), "resources", "web.geo")
        self.bad_config_file   = os.path.join(os.path.dirname(__file__), "resources", "bad.geo")
        self.output_path       = os.path.join(os.path.dirname(__file__), "resources")

    def test_no_config_file(self):
        with self.assertRaises(ValueError):
            ModflowOutput(self.resource, config_file=None, exe_name="mf2005", verbose=True)

    def test_bad_config_file(self):
        with self.assertRaises(ValueError):
            ModflowOutput(self.resource, config_file=self.bad_config_file, exe_name="mf2005", verbose=True)

    def test_load(self):
        mf = ModflowOutput(self.resource, config_file=self.wgs84_config_file, exe_name="mf2005", verbose=True)
        assert mf is not None

    def test_config_global_attributes(self):
        mf = ModflowOutput(self.resource, config_file=self.wgs84_config_file, exe_name="mf2005", verbose=True)
        assert mf is not None
        output_file = os.path.join(self.output_path, "convert.nc")
        mf.to_netcdf(output_file=output_file)
        nc = netCDF4.Dataset(output_file)
        self.assertEquals(nc.id, 'wgs84-config-test')
        self.assertEquals(nc.creator, 'modflow2netcdf')

    def test_plot(self):
        mf = ModflowOutput(self.resource, config_file=self.wgs84_config_file, exe_name="mf2005", verbose=True)
        assert mf is not None
        mf.to_plot().show()

    def test_web_config_file(self):
        mf = ModflowOutput(self.resource, config_file=self.web_config_file, exe_name="mf2005", verbose=True)
        assert mf is not None
        output_file = os.path.join(self.output_path, "web_config.nc")
        mf.to_netcdf(output_file=output_file)
        nc = netCDF4.Dataset(output_file)
        assert nc is not None

    def test_convert(self):
        mf = ModflowOutput(self.resource, config_file=self.wgs84_config_file, exe_name="mf2005", verbose=True)
        assert mf is not None
        output_file = os.path.join(self.output_path, "convert.nc")
        mf.to_netcdf(output_file=output_file)
        nc = netCDF4.Dataset(output_file)
        assert nc is not None

    def test_colorado_netcdf(self):
        working_dir = os.path.join(self.output_path, "colorado")
        colorado_nam = os.path.join(working_dir, "mod16_ssfix_wel2.nam")
        colorado_geo = os.path.join(working_dir, "colorado.geo")
        mf = ModflowOutput(colorado_nam, config_file=colorado_geo, exe_name="mf2005", verbose=True, model_ws=working_dir)
        assert mf is not None
        output_file = os.path.join(self.output_path, "colorado", "colorado.nc")
        mf.to_netcdf(output_file=output_file)
        nc = netCDF4.Dataset(output_file)
        assert nc is not None
        assert nc.variables.get("time").units == "minutes since 1970-01-01T00:00:00Z"

    def test_colorado_plot(self):
        working_dir = os.path.join(self.output_path, "colorado")
        colorado_nam = os.path.join(working_dir, "mod16_ssfix_wel2.nam")
        colorado_geo = os.path.join(working_dir, "colorado.geo")
        mf = ModflowOutput(colorado_nam, config_file=colorado_geo, exe_name="mf2005", verbose=True, model_ws=working_dir)
        assert mf is not None
        mf.to_plot(variable='heads', level=4).show()

    def test_freyberg_netcdf(self):
        working_dir = os.path.join(self.output_path, "freyberg")
        freyberg_nam = os.path.join(working_dir, "freyberg.nam")
        freyberg_geo = os.path.join(working_dir, "freyberg.geo")
        mf = ModflowOutput(freyberg_nam, config_file=freyberg_geo, exe_name="mf2005", verbose=True, model_ws=working_dir)
        assert mf is not None
        output_file = os.path.join(self.output_path, "freyberg", "freyberg.nc")
        mf.to_netcdf(output_file=output_file)
        nc = netCDF4.Dataset(output_file)
        assert nc is not None
        assert nc.variables.get("time").units == "days since 1970-01-01T00:00:00Z"

    def test_freyberg_plot(self):
        working_dir = os.path.join(self.output_path, "freyberg")
        freyberg_nam = os.path.join(working_dir, "freyberg.nam")
        freyberg_geo = os.path.join(working_dir, "freyberg.geo")
        mf = ModflowOutput(freyberg_nam, config_file=freyberg_geo, exe_name="mf2005", verbose=True, model_ws=working_dir)
        assert mf is not None
        mf.to_plot(variable='flow_right_face_centered', colormap=matplotlib.cm.GnBu).show()

    def test_carolina_netcdf(self):
        working_dir = os.path.join(self.output_path, "carolina")
        carolina_nam = os.path.join(working_dir, "1-23-07_NOASP.nam")
        carolina_geo = os.path.join(working_dir, "carolina.geo")
        mf = ModflowOutput(carolina_nam, config_file=carolina_geo, exe_name="mf2005", verbose=True, model_ws=working_dir)
        assert mf is not None
        output_file = os.path.join(self.output_path, "carolina", "carolina.nc")
        mf.to_netcdf(output_file=output_file)
        nc = netCDF4.Dataset(output_file)
        assert nc is not None
        assert nc.variables.get("time").units == "days since 1899-12-31T00:00:00Z"

    def test_carolina_plot(self):
        working_dir = os.path.join(self.output_path, "carolina")
        carolina_nam = os.path.join(working_dir, "1-23-07_NOASP.nam")
        carolina_geo = os.path.join(working_dir, "carolina.geo")
        mf = ModflowOutput(carolina_nam, config_file=carolina_geo, exe_name="mf2005", verbose=True, model_ws=working_dir)
        assert mf is not None
        mf.to_plot(variable='heads', colormap=matplotlib.cm.GnBu).show()

    def test_miami_dade_netcdf(self):
        working_dir = os.path.join(self.output_path, "miami-dade")
        miami_nam = os.path.join(working_dir, "UMD_fb.nam")
        miami_geo = os.path.join(working_dir, "miami-dade.geo")
        mf = ModflowOutput(miami_nam, config_file=miami_geo, exe_name="mf2005", verbose=True, model_ws=working_dir)
        assert mf is not None
        output_file = os.path.join(working_dir, "miami-dade.nc")
        mf.to_netcdf(output_file=output_file)
        nc = netCDF4.Dataset(output_file)
        assert nc is not None
        assert nc.variables.get("time").units == "days since 1996-01-01T00:00:00Z"

    def test_miami_dade_plot(self):
        working_dir = os.path.join(self.output_path, "miami-dade")
        miami_nam = os.path.join(working_dir, "UMD_fb.nam")
        miami_geo = os.path.join(working_dir, "miami-dade.geo")
        mf = ModflowOutput(miami_nam, config_file=miami_geo, exe_name="mf2005", verbose=True, model_ws=working_dir)
        assert mf is not None
        mf.to_plot(variable='heads', colormap=matplotlib.cm.GnBu).show()
