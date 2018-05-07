# -*- coding: utf-8 -*-

import numpy as np


from ...defaults.typing import *
from .exceptions import *


class GeoTransform(object):
    """
    Manage geotransform parameters for rasters.
    It is based on the GDAL GeoTransform concept.
    See: http://www.gdal.org/gdal_datamodel.html
    """

    # class constructor 
    def __init__(self, topLeftX: Number, topLeftY: Number, pixWidth: Number, pixHeight: Number, rotGT2: Number, rotGT4: Number) -> None:
        """
        Class constructor.

        :param topLeftX: top left corner of the top left pixel of the raster - x coord
        :type topLeftX: Number
        :param topLeftY:: top left corner of the top left pixel of the raster - y coord
        :type topLeftY: Number
        :param pixWidth: pixel width
        :type pixWidth: Number
        :param pixHeight: pixel height
        :type pixHeight: Number
        :param rotGT2: rotation
        :type rotGT2: Number
        :param rotGT4: rotation
        :type rotGT4: Number

        :return: None
        """

        self.gt = np.array([
            topLeftX,   # GT(0) - top left corner of the top left pixel of the raster
            pixWidth,   # GT(1) - pixel width
            rotGT2,     # GT(2) - rotation
            topLeftY,   # GT(3) - top left corner of the top left pixel of the raster
            rotGT4,     # GT(4) - rotation
            pixHeight   # GT(5) - pixel height
        ])

    @property
    def topLeftX(self):
        """
        Get top-left corner x value of the rasters.

        @return:  the top-left corner x value, according to GDAL convention - float.
        """

        return self.gt[0]

    @property
    def topLeftY(self):
        """
        Get top-left corner y value of the rasters.

        @return:  the top-left corner y value, according to GDAL convention - float.
        """

        return self.gt[3]

    @property
    def pixWidth(self):
        """
        Get East-West size of the rasters cell.

        @return:  the East-West size of the rasters cell - float.
        """

        return self.gt[1]

    @property
    def pixHeight(self):
        """
        Get North-South size of the rasters cell.

        @return:  the North-South size of the rasters cell - float.
        """

        return self.gt[5]

    @property
    def rotGT2(self):
        """
        Get rotation GT(2) (see GDAL documentation).

        @return:  the rasters rotation value GT(2). - float.
        """

        return self.gt[2]

    @property
    def rotGT4(self):
        """
        Get rotation GT(4) (see GDAL documentation).

        @return:  the rasters rotation value GT(4) - float.
        """

        return self.gt[4]

    def projectXY(self, xPixel: Number, yLine: Number) -> Tuple[float, float]:
        """
        Transforms from raster to geographic coordinates.

        See: http://www.gdal.org/gdal_datamodel.html
        "Note that the pixel/line coordinates in the above are from (0.0,0.0)
        at the top left corner of the top left pixel to (width_in_pixels,
        height_in_pixels) at the bottom right corner of the bottom right pixel.
        The pixel/line location of the center of the top left pixel would
        therefore be (0.5,0.5)."

        :param xPixel:
        :param yLine:
        :return: tuple storing geographic x - y pair
        """

        Xgeo = self.gt[0] + xPixel * self.gt[1] + yLine * self.gt[2]
        Ygeo = self.gt[3] + xPixel * self.gt[4] + yLine * self.gt[5]

        return Xgeo, Ygeo

    def check_params(self, tolerance=1e-06):
        """]
        Check absence of axis rotations or pixel size differences in the rasters band.

        @param  tolerance:  the maximum threshold for both pixel N-S and E-W difference, or axis rotations.
        @type  tolerance:  float.

        @return:  None when successful, RasterParametersException when pixel differences or axis rotations.

        @raise: RasterParametersException - rasters geometry incompatible with this module (i.e. different cell sizes or axis rotations).
        """
        
        # check if pixel size can be considered the same in the two axis directions
        if abs(abs(self._pixWidth) - abs(self._pixHeight)) / abs(self._pixHeight) > tolerance:
            raise RasterParametersException("Pixel sizes in x and y directions are different in rasters")

            # check for the absence of axis rotations
        if abs(self._rotation_GT_2) > tolerance or abs(self._rotation_GT_4) > tolerance:
            raise RasterParametersException("There should be no axis rotation in rasters")

        return

    def llcorner(self):
        """
        Creates a point tuple (x, y) representing the lower-left corner of the rasters.

        @return:  new Point instance.
        """
        return self.topLeftX, self.topLeftY - abs(self.pixSizeNS) * self.rows

    def trcorner(self):
        """
        Create a point tuple (x, y) representing the top-right corner of the rasters.

        @return:  new Point instance.
        """
        return self.topLeftX + abs(self.pixSizeEW) * self.cols, self.topLeftY

    def geo_equiv(self, other, tolerance=1.0e-6):
        """
        Checks if two rasters are geographically equivalent.

        @param  other:  a grid to be compared with self.
        @type  other:  Raster instance.
        @param  tolerance:  the maximum threshold for pixel sizes, topLeftX or topLeftY differences.
        @type  tolerance:  float.

        @return:  Boolean.
        """
        if 2 * (self.topLeftX - other.topLeftX) / (self.topLeftX + other.topLeftX) > tolerance or \
                                        2 * (self.topLeftY - other.topLeftY) / (
                                    self.topLeftY + other.topLeftY) > tolerance or \
                                        2 * (abs(self.pixSizeEW) - abs(other.pixSizeEW)) / (
                                    abs(self.pixSizeEW) + abs(other.pixSizeEW)) > tolerance or \
                                        2 * (abs(self.pixSizeNS) - abs(other.pixSizeNS)) / (
                                    abs(self.pixSizeNS) + abs(other.pixSizeNS)) > tolerance or \
                        self.rows != other.rows or self.cols != other.cols or self.projection != other.projection:
            return False
        else:
            return True

