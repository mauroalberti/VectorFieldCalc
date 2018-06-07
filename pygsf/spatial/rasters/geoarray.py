# -*- coding: utf-8 -*-


from ...mathematics.defaults import *
from ...defaults.typing import *

from .exceptions import *

from ...mathematics.arrays import interp_bilinear
from pygsf.spatial.rasters.geotransform import GeoTransform, geogrToPix, pixToGeogr

from .fields import orients_d, magnitude, magn_grads, divergence, curl_module, magn_grad_along_flowlines


class GeoArray(object):
    """
    GeoArray class.
    Stores and process georeferenced raster data.

    """

    def __init__(self, inGeotransform: GeoTransform, inProjection: str, inLevels: Optional[List['np.array']]=None) -> None:
        """
        GeoArray class constructor.

        :param  inGeotransform:  the geotransform
        :type  inGeotransform:  GeoTransform.
        :param inProjection: the projection
        :type inProjection: str
        :param  inLevels:  the nd-array storing the data.
        :type  inLevels:  np.array.

        :return:  None.

        Examples:
        """

        self.gt = inGeotransform
        self.prj = inProjection
        if inLevels is None:
            self._levels = []
        else:
            self._levels = inLevels

    @property
    def cellsize_x(self) -> float:
        """
        Get the cell size of the grid in the x direction.

        :return: cell size in the x (j) direction.
        :rtype: float.

        Examples:
        """

        return abs(self.gt.pixWidth)

    @property
    def cellsize_y(self) -> float:
        """
        Get the cell size of the grid in the y direction.

        :return: cell size in the y (-i) direction.
        :rtype: float.

        Examples:
        """

        return abs(self.gt.pixHeight)

    @property
    def levels_num(self) -> int:
        """
        Returns the number of levels (dimensions) of the grid.

        :return: number of levels.
        :rtype: int.

        Examples:
          >>> GeoArray(grid_data=array([[1, 2], [3, 4]])).levels_num
          1
          >>> GeoArray(grid_data=np.ones((4, 3, 2)))
          2
        """

        return len(self._levels)

    def level(self, level_ndx: int=0) -> Optional['array']:
        """
        Return the array corresponding to the requested level
        if existing else None.

        :param level_ndx: the index of the requested level.
        :type level_ndx: int.
        :return: the array or None.
        :rtype: optional array.

        Examples:
        """

        if 0 <= level_ndx < self.levels_num:
            return self._levels[level_ndx]
        else:
            return None

    def level_shape(self, level_ndx: int=0) -> Optional[Tuple[int, int]]:
        """
        Returns the shape (num. rows and num. columns) of the considered level grid.

        :param level_ndx: index of the level (grid) to consider.
        :type level_ndx: int.
        :return: number of rows and columns of the specific grid.
        :rtype: optional tuple of two int values.

        Examples:
          >>> GeoArray(grid_data=array([[1, 2], [3, 4]])).level_shape
          2, 2
          >>> GeoArray(grid_data=np.ones((4, 3, 2)))
          1
        """

        if 0 <= level_ndx < self.levels_num:
            return self._levels[level_ndx].shape
        else:
            return None

    def level_llc(self, level_ndx: int=0) -> Optional[Tuple[int, int]]:
        """
        Returns the coordinates of the lower-left corner.

        :param level_ndx: index of the level (grid) to consider.
        :type level_ndx: int.
        :return: x and y values of the lower-left corner of the specific grid.
        :rtype: optional tuple of two int values.

        Examples:
        """

        shape = self.level_shape(level_ndx)
        if not shape:
            return None

        llc_i, llc_j = shape[0], 0

        return self.ijToxy(llc_i, llc_j)

    def xyToij(self, x: Number, y: Number) -> Tuple[Number, Number]:
        """
        Converts from geographic to array coordinates.

        :param x: x geographic component
        :type x: Number
        :param y: y geographic component
        :type y: Number
        :return: i and j values
        :type: tuple of two float values

        Examples:
        """

        return geogrToPix(self.gt, x, y)

    def ijToxy(self, i: Number, j: Number) -> Tuple[Number, Number]:
        """
        Converts from array to geographic coordinates.

        :param i: i array component.
        :type i: Number.
        :param j: j array component.
        :type j: Number.
        :return: x and y geographic coordinates.
        :type: tuple of two float values.

        Examples:
        """

        return pixToGeogr(self.gt, i, j)

    @property
    def has_rotation(self) -> bool:
        """
        Determines if a geoarray has axis rotations defined.

        :return: true if there are rotations, false otherwise.
        :rtype: bool.

        Examples:
        """

        return self.gt.rotRow != 0.0 or self.gt.rotColumn != 0.0

    def interpolate_bilinear(self, x: Number, y: Number, level_ndx=0) -> Optional[float]:
        """
        Interpolate the z value at a point, given its geographic coordinates.
        Interpolation method: bilinear.

        :param x: x geographic coordinate.
        :type x: Number.
        :param y: y geographic coordinate.
        :type y: Number.
        :param level_ndx: the index of the used array.
        :type level_ndx: int.
        :return: a geoarray storing the interpolated z value
        :rtype: optional float.

        Examples:
        """

        i, j = self.xyToij(x, y)

        return interp_bilinear(self._levels[level_ndx], i, j)

    def magnitude_field(self, ndx_fx=0, ndx_fy=1) -> 'GeoArray':
        """
        Calculates magnitude field as a geoarray.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the magnitude field.
        :rtype: GeoArray.

        Examples:
        """

        magn = magnitude(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy])

        return GeoArray(
            inGeotransform=self.gt,
            inProjection=self.prj,
            inLevels=[magn]
        )

    def orientations(self, ndx_fx=0, ndx_fy=1) -> 'GeoArray':
        """
        Calculates orientations field as a geoarray.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the orientation field.
        :rtype: GeoArray.

        Examples:
        """

        orient = orients_d(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy])

        return GeoArray(
            inGeotransform=self.gt,
            inProjection=self.prj,
            inLevels=[orient]
        )

    def divergence_2D(self, ndx_fx=0, ndx_fy=1) -> 'GeoArray':
        """
        Calculates divergence of a 2D field as a geoarray.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the divergence field.
        :rtype: GeoArray.

        Examples:
        """

        div = divergence(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy],
            cell_size_x=self.cellsize_x,
            cell_size_y=self.cellsize_y)

        return GeoArray(
            inGeotransform=self.gt,
            inProjection=self.prj,
            inLevels=[div]
        )

    def curl_module(self, ndx_fx=0, ndx_fy=1) -> 'GeoArray':
        """
        Calculates curl module of a 2D field as a geoarray.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the curl module field.
        :rtype: GeoArray.

        Examples:
        """

        curl_m = curl_module(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy],
            cell_size_x=self.cellsize_x,
            cell_size_y=self.cellsize_y)

        return GeoArray(
            inGeotransform=self.gt,
            inProjection=self.prj,
            inLevels=[curl_m])

    def magnitude_grads(self, axis: str= '', ndx_fx: int=0, ndx_fy: int=1) -> 'GeoArray':
        """
        Calculates the magnitude gradient along the x axis of a 2D field as a geoarray.

        :param axis: axis along wich to calculate the gradient, 'x' or 'y', or '' (predefined) for both x and y.
        :type axis: str.
        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the magnitude gradient along the x axis field.
        :rtype: GeoArray.
        :raises: GeoArrayIOException.

        Examples:
        """

        if axis == 'x':
            cell_sizes = [self.cellsize_x]
        elif axis == 'y':
            cell_sizes = [self.cellsize_y]
        elif axis == '':
            cell_sizes = [self.cellsize_x, self.cellsize_y]
        else:
            raise GeoArrayIOException("Axis must be 'x' or 'y. '{}' given".format(axis))

        magn_grads = magn_grads(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy],
            dir_cell_sizes=cell_sizes,
            axis=axis)

        return GeoArray(
            inGeotransform=self.gt,
            inProjection=self.prj,
            inLevels=magn_grads)

    def grad_flowlines(self, ndx_fx: int=0, ndx_fy: int=1) -> 'GeoArray':
        """
        Calculates gradient along flow lines.

        :param ndx_fx: index of x field.
        :type ndx_fx: integer.
        :param ndx_fy: index of y field.
        :type ndx_fy: integer.
        :return: a geoarray storing the flowline gradient field
        :rtype: GeoArray
        """

        flowln_grad = magn_grad_along_flowlines(
            fld_x=self._levels[ndx_fx],
            fld_y=self._levels[ndx_fy],
            cell_size_x=self.cellsize_x,
            cell_size_y=self.cellsize_y)

        return GeoArray(
            inGeotransform=self.gt,
            inProjection=self.prj,
            inLevels=[flowln_grad])


if __name__ == "__main__":

    import doctest
    doctest.testmod()
