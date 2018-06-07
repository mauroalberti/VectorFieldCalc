# coding: utf-8

# # Divergence result check

# Theoretical example based on:
# 
# https://www.khanacademy.org/math/multivariable-calculus/multivariable-derivatives/divergence-and-curl-articles/a/divergence
# 

# ## Required imports

import unittest

from ..spatial.rasters.fields import *


# Z transfer functions

# These analytical functions define the value of the cells, from the provided x and y geographic coordinates.

def z_func_fx(x, y):
    return 2 * x - y


def z_func_fy(x, y):
    return y * y


def z_func_div(x, y):
    return 2 + 2 * y


class TestDivergence(unittest.TestCase):

    def setUp(self):

        pass

    def test_divergence(self):
        """
        Test the divergence calculation.

        :return:
        """

        # The geotransform defines the raster-to-geographic coordinates mapping.

        gt1 = GeoTransform(1500, 3000, 10, 10)

        # Z fields (x and y components)

        # The gridded field values are calculated for the vector field x- and y- components,
        # as well as the expected teorethical divergence field.

        rows = 100
        cols = 50

        # Vector field x-component

        fx = array_from_function(row_num=rows, col_num=cols, geotransform=gt1, z_transfer_func=z_func_fx)

        # Vector field y-component

        fy = array_from_function(row_num=rows, col_num=cols, geotransform=gt1, z_transfer_func=z_func_fy)

        # Theoretical divergence

        div_theor = array_from_function(row_num=rows, col_num=cols, geotransform=gt1, z_transfer_func=z_func_div)

        # Divergence as resulting from pygsf calculation:

        div_pygsf = divergence(fx, fy, 10, 10)

        assert np.allclose(div_theor, div_pygsf)

