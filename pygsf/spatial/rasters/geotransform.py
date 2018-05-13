# -*- coding: utf-8 -*-


import numpy as np

from ...defaults.typing import *


class GeoTransform(object):
    """
    Manage geotransform parameters for rasters.
    It is based on the GDAL GeoTransform concept.
    See: http://www.gdal.org/gdal_datamodel.html
    """

    def __init__(self, inTopLeftX: Number, inTopLeftY: Number, inPixWidth: Number, inPixHeight: Number, inRotRow: Number=0.0, inRotColumn: Number=0.0) -> None:
        """
        Class constructor.

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

        self.gt = np.array([
            inTopLeftX,   # GT(0) - top left corner of the top left pixel of the raster
            inPixWidth,   # GT(1) - pixel width
            inRotRow,     # GT(2) - row rotation
            inTopLeftY,   # GT(3) - top left corner of the top left pixel of the raster
            inRotColumn,  # GT(4) - column rotation
            inPixHeight   # GT(5) - pixel height
        ], dtype=float)

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

        return self.gt[0], self.gt[1], self.gt[2], self.gt[3], self.gt[4], self.gt[5]

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

        return self.gt[0]

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

        return self.gt[3]

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

        return self.gt[1]

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

        return self.gt[5]

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

        return self.gt[2]

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

        return self.gt[4]

    def pixToGeogr(self, xPixel: Number, yLine: Number) -> Tuple[float, float]:
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
          >>> gt1.pixToGeogr(10, 10)
          (1600.0, 3100.0)
        """

        Xgeo = self.gt[0] + xPixel * self.gt[1] + yLine * self.gt[2]
        Ygeo = self.gt[3] + xPixel * self.gt[4] + yLine * self.gt[5]

        return Xgeo, Ygeo

    def geogrToPix(self, x: Number, y: Number) -> Tuple[float, float]:
        """
        Transforms from geographic to pixel coordinates.

        See: http://www.gdal.org/gdal_datamodel.html
        "Note that the pixel/line coordinates in the above are from (0.0,0.0)
        at the top left corner of the top left pixel to (width_in_pixels,
        height_in_pixels) at the bottom right corner of the bottom right pixel.
        The pixel/line location of the center of the top left pixel would
        therefore be (0.5,0.5)."

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
          >>> gt1.geogrToPix(1600, 3100)
          (10.0, 10.0)
        """

        g0, g1, g2, g3, g4, g5 = self.components

        l = (g1*g3 - g0*g4 + g4*x - g1*y) / (g2*g4 - g1*g5)
        p = (x - g0 - l*g2) / g1

        return p, l


if __name__ == "__main__":

    import doctest
    doctest.testmod()

