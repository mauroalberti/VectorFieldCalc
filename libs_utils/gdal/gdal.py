# -*- coding: utf-8 -*-


import numpy as np

from typing import Tuple

import gdal
from gdalconst import *

from .exceptions import *

from pygsf.spatial.rasters.geotransform import GeoTransform


def read_raster(raster_path: str):
    """
    Read a raster layer.

    :param raster_path:
    :return:

    Examples:
    """

    # open raster file and check operation success

    raster_data = gdal.Open(str(raster_path), GA_ReadOnly)
    if raster_data is None:
        raise RasterIOException("No input data open")

    # set grid values from geotransform array

    raster_params.set_topLeftX(raster_data.GetGeoTransform()[0])
    raster_params.set_pixSizeEW(raster_data.GetGeoTransform()[1])
    raster_params.set_rotationA(raster_data.GetGeoTransform()[2])
    raster_params.set_topLeftY(raster_data.GetGeoTransform()[3])
    raster_params.set_rotationB(raster_data.GetGeoTransform()[4])
    raster_params.set_pixSizeNS(raster_data.GetGeoTransform()[5])



def read_raster_band(raster_name: str, bnd_ndx: int=1):
    """
    Read data and metadata of a rasters band based on GDAL.

    :param raster_name:
    :param raster_params:
    :return:
    """

    # get single band

    band = raster_data.GetRasterBand(1)

    # get no data value for current band

    raster_params.set_noDataValue(band.GetNoDataValue())
    if raster_params.get_noDataValue() is None:
        raise RasterIOException("Unable to get no data value from input rasters. Try change input format\n(e.g., ESRI ascii grids generally work)")

    # read data from band

    grid_values = band.ReadAsArray(0, 0, raster_params.get_cols(), raster_params.get_rows())
    if grid_values is None:
        raise RasterIOException("Unable to read data from rasters")

    # transform data into numpy array

    data = np.asarray(grid_values)

    # if nodatavalue exists, set null values to NaN in numpy array

    if raster_params.get_noDataValue() is not None:
        data = np.where(abs(data - raster_params.get_noDataValue()) > 1e-05, data, np.NaN)

    return raster_params, data


def read_raster_layer(raster_name: str, layermap_items) -> Tuple[GeoTransform, 'np.array']:
    """

    :param raster_name:
    :param layermap_items:
    :return:
    """

    # verify input parameters
    if raster_name is None or raster_name == '':
        raise RasterIOException("No name defined for rasters")

    # get rasters input file
    raster_layer = None
    for (name, layer) in layermap_items:
        if layer.name() == raster_name:
            raster_layer = layer
            break
    if raster_layer is None:
        raise RasterIOException("Unable to get rasters name")

    try:
        raster_source = raster_layer.source()
    except:
        raise RasterIOException("Unable to get rasters file")

    # get rasters parameters and data
    try:
        raster_params, raster_array = read_raster_band(raster_source)
    except Exception as e:
        raise RasterIOException(str(e))

    return raster_params, raster_array

