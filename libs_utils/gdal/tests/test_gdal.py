# -*- coding: utf-8 -*-


import unittest

from ..gdal import *


data_path = r"./example_data/vx.asc"


class TestRasterRead(unittest.TestCase):

    def setUp(self):

        pass

    def test_raster_read(self):

        dataset, geotransform, num_bands, projection = read_raster(data_path)

        assert num_bands == 1
        assert projection == ''

        dataset = None

    def test_band_read(self):

        dataset, _, _, _ = read_raster(data_path)

        band = read_band(dataset, 1)

        dataset = None

if __name__ == '__main__':

    unittest.main()

