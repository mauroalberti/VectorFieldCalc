# -*- coding: utf-8 -*-


import os

from .geoarray import GeoArray
from .exceptions import RasterIOExceptions


def write_esrigrid(geoarray: GeoArray, outgrid_fn: str, esri_nullvalue: Number=-99999, level: int=0):
    """
    writes ESRI ascii grid
    
    :param geoarray: 
    :param outgrid_fn: 
    :param esri_nullvalue: 
    :param level: 
    :return: 
    """
    
    outgrid_fn = str(outgrid_fn)

    # checking existence of output slope grid
    if os.path.exists(outgrid_fn):
        raise RasterIOExceptions("Output grid '%s' already exists" % outgrid_fn)

    try:
        outputgrid = open(outgrid_fn, 'w')  # create the output ascii file
    except:
        raise RasterIOExceptions("Unable to create output grid '%s'" % outgrid_fn)

    if outputgrid is None:
        raise RasterIOExceptions("Unable to create output grid '%s'" % outgrid_fn)

    # writes header of grid ascii file
    outputgrid.write("NCOLS %d\n" % geoarray.get_xlines_num())
    outputgrid.write("NROWS %d\n" % geoarray.get_ylines_num())
    outputgrid.write("XLLCORNER %.8f\n" % geoarray.grid_domain.get_start_point().x)
    outputgrid.write("YLLCORNER %.8f\n" % geoarray.grid_domain.get_start_point().y)
    outputgrid.write("CELLSIZE %.8f\n" % geoarray.get_cellsize_horiz_mean())
    outputgrid.write("NODATA_VALUE %d\n" % esri_nullvalue)

    esrigrid_outvalues = np.where(np.isnan(geoarray.grid_data), esri_nullvalue, geoarray.grid_data)

    # output of results
    if len(geoarray.grid_data.shape) == 3:
        for i in range(0, geoarray.get_ylines_num()):
            for j in range(0, geoarray.get_xlines_num()):
                outputgrid.write("%.8f " % (esrigrid_outvalues[i, j, level]))
            outputgrid.write("\n")
    elif len(geoarray.grid_data.shape) == 2:
        for i in range(0, geoarray.get_ylines_num()):
            for j in range(0, geoarray.get_xlines_num()):
                outputgrid.write("%.8f " % (esrigrid_outvalues[i, j]))
            outputgrid.write("\n")

    outputgrid.close()

    return True

