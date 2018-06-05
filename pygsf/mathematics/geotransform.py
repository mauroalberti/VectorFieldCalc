# -*- coding: utf-8 -*-


import numpy as np

from ..defaults.typing import *


class GeoTransform(np.ndarray):
    """
    Manage geotransform parameters for rasters.
    It is based on the GDAL GeoTransform concept.
    See: http://www.gdal.org/gdal_datamodel.html
    """

    def __new__(cls, inTopLeftX: Number, inTopLeftY: Number, inPixWidth: Number, inPixHeight: Number, inRotRow: Number=0.0, inRotColumn: Number=0.0) -> None:
        """
        Instance creator.
        Note: pixel height input is positive.

        :param inTopLeftX: top left corner of the top left pixel of the raster - x coord
        :type inTopLeftX: Number
        :param inTopLeftY:: top left corner of the top left pixel of the raster - y coord
        :type inTopLeftY: Number
        :param inPixWidth: pixel width
        :type inPixWidth: Number
        :param inPixHeight: pixel height
        :type inPixHeight: Number
        :param inRotRow: rotation
        :type inRotRow: Number
        :param inRotColumn: rotation
        :type inRotColumn: Number

        :return: None

        Examples:
          >>> GeoTransform(1500, 3000, 10, 10)
          GeoTransform(topLeftX: 1500.00, topLeftY: 3000.00, pixWidth: 10.00, pixHeight: 10.00, rotRow: 0.00, rotColumn: 0.00)
          """

        return np.array([
            inTopLeftX,  # GT(0) - top left corner of the top left pixel of the raster
            inPixWidth,  # GT(1) - pixel width
            inRotRow,  # GT(2) - row rotation
            inTopLeftY,  # GT(3) - top left corner of the top left pixel of the raster
            inRotColumn,  # GT(4) - column rotation
            -inPixHeight  # GT(5) - pixel height
        ], dtype=float).view(cls)

    @classmethod
    def fromGdalGt(cls, gdal_gt: Tuple[float, float, float, float, float, float]):
        """
        Creates a Geotransform from a GDAL-convetion tuple.

        :param gdal_gt: tuple with GDAL geotransform parameters
        :return: None
        """

        return cls(
            gdal_gt[0],
            gdal_gt[3],
            gdal_gt[1],
            gdal_gt[5],
            gdal_gt[2],
            gdal_gt[4])

    def __repr__(self) -> str:

        return "GeoTransform(topLeftX: {:.2f}, topLeftY: {:.2f}, pixWidth: {:.2f}, pixHeight: {:.2f}, rotRow: {:.2f}, rotColumn: {:.2f})".format(
            self.topLeftX,
            self.topLeftY,
            self.pixWidth,
            self.pixHeight,
            self.rotRow,
            self.rotColumn)

    @property
    def components(self):
        """
        Returns the Geotransform components as a tuple.

        :return:
        """

        return self[0], self[1], self[2], self[3], self[4], self[5]

    @property
    def topLeftX(self) -> float:
        """
        Get top-left corner x value of the rasters.

        :return: the top-left corner x value, according to GDAL convention
        :rtype: float.

        Examples:
          >>> GeoTransform(1500, 3000, 10, 10, 0, 0).topLeftX
          1500.0
        """

        return self[0]

    @property
    def topLeftY(self) -> float:
        """
        Get top-left corner y value of the rasters.

        :return:  the top-left corner y value, according to GDAL convention.
        :rtype: float.

        Examples:
          >>> GeoTransform(1500, 3000, 10, 10, 0, 0).topLeftY
          3000.0
          """

        return self[3]

    @property
    def pixWidth(self) -> float:
        """
        Get East-West size of the rasters cell.

        :return:  the East-West size of the rasters cell
        :rtype: float.

        Examples:
          >>> GeoTransform(1500, 3000, 10, 10, 0, 0).pixWidth
          10.0
        """

        return self[1]

    @property
    def pixHeight(self) -> float:
        """
        Get North-South size of the rasters cell.

        :return:  the North-South size of the rasters cell.
        :rtype: float.

        Examples:
          >>> GeoTransform(1500, 3000, 10, 10, 0, 0).pixHeight
          10.0
        """

        return -self[5]

    @property
    def rotRow(self) -> float:
        """
        Get row rotation GT(2) (see GDAL documentation).

        :return:  the rasters rotation value GT(2).
        :rtype: float.

        Examples:
          >>> GeoTransform(1500, 3000, 10, 10, 0, 0).rotRow
          0.0
        """

        return self[2]

    @property
    def rotColumn(self) -> float:
        """
        Get column rotation GT(4) (see GDAL documentation).

        :return:  the rasters rotation value GT(4).
        :rtype: float.

        Examples:
          >>> GeoTransform(1500, 3000, 10, 10, 0, 0).rotColumn
          0.0
        """

        return self[4]


def pixToGeogr(geotransform: GeoTransform, i: Number, j: Number) -> Tuple[float, float]:
    """
    Transforms from pixel to geographic coordinates.

    See: http://www.gdal.org/gdal_datamodel.html
    "Note that the pixel/line coordinates in the above are from (0.0,0.0)
    at the top left corner of the top left pixel to (width_in_pixels,
    height_in_pixels) at the bottom right corner of the bottom right pixel.
    The pixel/line location of the center of the top left pixel would
    therefore be (0.5,0.5)."

    :param xPixel: the  pixel x coordinate.
    :type xPixel: Number.
    :param yLine: the pixel y coordinate.
    :type yLine: Number
    :return: tuple storing geographic x-y pair
    :rtype: tuple of two floats.

    Examples:
      >>> gt1 = GeoTransform(1500, 3000, 10, 10)
      >>> pixToGeogr(gt1, 0, 0)
      (1500.0, 3000.0)
      >>> pixToGeogr(gt1, 10, 10)
      (1600.0, 2900.0)
      >>> pixToGeogr(gt1, 2, 1)
      (1510.0, 2980.0)
    """

    Xgeo = geotransform[0] + j * geotransform[1] + i * geotransform[2]
    Ygeo = geotransform[3] + j * geotransform[4] + i * geotransform[5]

    return Xgeo, Ygeo


def geogrToPix(geotransform: GeoTransform, x: Number, y: Number) -> Tuple[float, float]:
    """
    Transforms from geographic to pixel coordinates.

    See: http://www.gdal.org/gdal_datamodel.html
    "Note that the pixel/line coordinates in the above are from (0.0,0.0)
    at the top left corner of the top left pixel to (width_in_pixels,
    height_in_pixels) at the bottom right corner of the bottom right pixel.
    The pixel/line location of the center of the top left pixel would
    therefore be (0.5,0.5)."

    Formula derivation
    ------------------

    x = g0 + p * g1 + l * g2
    y = g3 + p * g4 + l * g5

    x - g0 - l* g2 = p * g1
    p = (x - g0 - l*g2) / g1
    p = (y - g3 - l*g5) / g4

    g1 * y - g1 * g3 - l * g5 * g1 = g4 * x - g0 * g4 - l * g2 * g4
    l * (g2g4 -g1g5) = -g1y + g1g3 + g4x - g0g4
    l = (g1g3 - g0g4 + g4x - g1y) / (g2g4 - g1g5)

    (x - g0 - g1p) / g2 =  (g1g3 - g0g4 + g4x - g1y) / (g2g4 - g1g5)

    g2 * (g1g3 - g0g4 + g4x - g1y) / (g2g4 - g1g5) = x - g0 - g1p


    :param x: the  geographic x coordinate.
    :type x: Number.
    :param y: the geographic y coordinate.
    :type y: Number
    :return: tuple storing pixel x-y pair
    :rtype: tuple of two floats.

    Examples:
      >>> gt1 = GeoTransform(1500, 3000, 10, 10)
      >>> geogrToPix(gt1, 1600, 2900)
      (10.0, 10.0)
    """

    g0, g1, g2, g3, g4, g5 = geotransform.components

    l = (g1*g3 - g0*g4 + g4*x - g1*y) / (g2*g4 - g1*g5)
    p = (x - g0 - l*g2) / g1

    return p, l


def gtToXY(gt: GeoTransform, num_rows: int, num_cols: int) -> Tuple['array', 'array']:
    """
    Create two arrays that represent the X and Y geographic coordinates of
    the cells CENTERS (not corners) given the geotransform.

    :param gt: the source geotransform.
    :type gt: GeoTransform
    :param num_rows: the number of rows.
    :type num_rows: int.
    :param num_cols: the number of the columns.
    :type num_cols: int.
    :return: the two Numpy arrays representing the geographic X and Y coordinates.
    :rtype: Numpy array of float64.

    Examples:
    """

    X = np.zeros((num_rows, num_cols), dtype=np.float64)
    Y = np.zeros((num_rows, num_cols), dtype=np.float64)
    for i in range(num_rows):
        for j in range(num_cols):
            x, y = pixToGeogr(gt, i+0.5, j+0.5)
            X[i, j] = x
            Y[i, j] = y

    return X, Y


if __name__ == "__main__":

    import doctest
    doctest.testmod()

