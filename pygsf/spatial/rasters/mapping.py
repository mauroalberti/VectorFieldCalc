# -*- coding: utf-8 -*-

from typing import List, Callable

import numpy as np

from .geotransform import GeoTransform, pixToGeogr
from ...defaults.typing import *


def ij_transfer_func(i: Number, j: Number, geotransform: GeoTransform, z_transfer_func: Callable, i_shift=0.5, j_shift=0.5):
    """
    Return a z value as the result of a function (transfer_func_z) applied to a
    (i+i_shift,j+j_shift) point (i.e., with defaultvalues, the cell center, not the cell top-left corner)
    given a geotransform.

    :param  i:  array i (-y) coordinate of a single point.
    :type  i:  Number.
    :param  j:  array j (x) coordinate of a single point.
    :type  j:  Number.
    :param  geotransform:  geotransform
    :type  geotransform:  GeoTransform.
    :param  z_transfer_func:  function that calculates the z value given x and y input
    :type  z_transfer_func:  function.
    :param i_shift: cell unit shift in the i direction with respect to the cell top-left corner. Default is 0.5, i.e. half the cell size
    :type i_shift: Number
    :param j_shift: cell unit shift in the j direction with respect to the cell top-left corner. Default is 0.5, i.e. half the cell size
    :type j_shift: Number
    :return: z value
    :rtype: float.

    Examples:
    """

    return z_transfer_func(*pixToGeogr(geotransform, i+i_shift, j+j_shift))


def array_from_function(row_num: int, col_num: int, geotransform: GeoTransform, z_transfer_func: Callable) -> 'np.array':
    """
    Creates an array of z values based on functions that map (i,j) indices (to be created)
    into (x, y) values and then z values.

    :param  row_num:  row number of the array to be created.
    :type  row_num:  int.
    :param  col_num:  column number of the array to be created.
    :type  col_num:  int.
    :param  geotransform:  geotransform.
    :type  geotransform:  GeoTransform.
    :param  z_transfer_func:  function that derives z given a (x, y) point.
    :type  z_transfer_func:  Function.

    :return:  array of z values
    :rtype: np.array of float numbers.

    Examples:
    """

    return np.fromfunction(
        function=ij_transfer_func,
        shape=(row_num, col_num),
        geotransform=geotransform,
        z_transfer_func=z_transfer_func)


if __name__ == "__main__":

    import doctest
    doctest.testmod()
