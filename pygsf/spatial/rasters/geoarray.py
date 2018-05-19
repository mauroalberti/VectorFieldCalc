# -*- coding: utf-8 -*-


import copy

from ...mathematics.defaults import *
from ...defaults.typing import *

from .exceptions import *

from .geotransform import GeoTransform


class GeoArray(object):
    """
    GeoArray class.
    Stores and process georeferenced raster data.

    """

    def __init__(self, inGeotransform: GeoTransform, inProjection: str, inArray: 'array') -> None:
        """
        GeoArray class constructor.

        :param  inGeotransform:  the geotransform
        :type  inGeotransform:  GeoTransform.
        :param inProjection: the projection
        :type inProjection: str
        :param  inArray:  the nd-array storing the data.
        :type  inArray:  np.array.

        :return:  None.

        Examples:
        """

        self.gt = inGeotransform
        self.prj = inProjection
        self.arr = inArray

    @property
    def row_num(self) -> int:
        """
        Get row number of the grid domain.

        :return: number of rows of data array.
        :rtype: int.

        Examples:
        """

        return np.shape(self.arr)[0]

    @property
    def col_num(self) -> int:
        """
        Get column number of the grid domain.

        :return: number of columns of data array
        :rtype: int.

        Examples:
        """

        return np.shape(self.arr)[1]

    @property
    def cellsize_x(self) -> float:
        """
        Get the cell size of the grid in the x direction.

        :return: cell size in the x (j) direction
        :rtype: float.

        Examples:
        """

        return abs(self.gt.pixWidth)

    @property
    def cellsize_y(self) -> float:
        """
        Get the cell size of the grid in the y direction.

        :return: cell size in the y (-i) direction
        :rtype: float.

        Examples:
        """

        return abs(self.gt.pixHeight)

    @property
    def cellsize_h(self) -> float:
        """
        Get the mean horizontal cell size.

        :return: mean horizontal cell size
        :rtype: float.

        Examples:
        """

        return (self.cellsize_x + self.cellsize_y) / 2.0

    def geog2array_coord(self, x: float, y: float) -> Tuple[float, float]:
        """
        Converts from geographic to rasters (array) coordinates.

        :param x: x coordinate of point to sample.
        :type x: float.
        :param y: ycoordinate of point to sample.
        :type y: float.
        :return: point coordinates in rasters (array) frame
        :rtype: ArrCoord.

        Examples:
        """

        return self.gt.geogrToPix(x, y)

    def x(self) -> 'array':
        """
        Creates an array storing the geographical coordinates of the cell centers along the x axis.
        Direction is from left to right.

        :return: array with shape: 1 x col_num.
        :rtype: numpy.array

        Examples:
        """

        x_values = self.domain.llcorner.x + self.cellsize_x * (0.5 + np.arange(self.col_num))

        return x_values[np.newaxis, :]

    def y(self) -> 'array':
        """
        Creates an array storing the geographical coordinates of the cell centers along the y axis.
        Direction is from top to bottom.

        :return: array with shape: 1 x col_num.
        :rtype: numpy.array

        Examples:
        """

        y_values = self.domain.trcorner.y - self.cellsize_y * (0.5 + np.arange(self.row_num))

        return y_values[:, np.newaxis]

    @property
    def levels_num(self) -> int:
        """
        Returns the number of levels (dimensions) of the grid.

        :return: number of levels
        :rtype: int

        Examples:
          >>> GeoArray(grid_data=array([[1, 2], [3, 4]])).levels_num
          1
          >>> GeoArray(grid_data=np.ones((4, 3, 2)))
          2
        """

        return self.arr.ndim

    def grad_forward_y(self) -> 'array':
        """
        Return an array representing the forward gradient in the y direction (top-wards), with values scaled by cell size.

        :return: array with same shape as current GeoArray instance
        :rtype: numpy.array

        Examples:
        """

        gf = np.zeros(np.shape(self.arr)) * np.NaN
        gf[1:, :] = self.arr[:-1, :] - self.arr[1:, :]

        return gf / float(self.cellsize_y)

    def grad_forward_x(self) -> 'array':
        """
        Return an array representing the forward gradient in the x direction (right-wards), with values scaled by cell size.

        :return: array with same shape as current GeoArray instance
        :rtype: numpy.array

        Examples:
        """

        gf = np.zeros(np.shape(self.arr), ) * np.NaN
        gf[:, :-1] = self.arr[:, 1:] - self.arr[:, :-1]

        return gf / float(self.cellsize_x)

    def interpolate_bilinear(self, array_coords: ArrCoord) -> float:
        """
        Interpolate the z value at a point, given its array coordinates.
        Interpolation method: bilinear.

        :param array_coords: array coordinates of the point for which the interpolation will be made.
        :type array_coords: ArrCoord.

        :return: interpolated z value
        :rtype: float.
        """

        loc_cellcenter_i = array_coords.i - 0.5
        loc_cellcenter_j = array_coords.j - 0.5

        assert loc_cellcenter_i > 0, loc_cellcenter_j > 0

        grid_val_00 = self.arr[int(floor(loc_cellcenter_i)), int(floor(loc_cellcenter_j))]
        grid_val_01 = self.arr[int(floor(loc_cellcenter_i)), int(ceil(loc_cellcenter_j))]
        grid_val_10 = self.arr[int(ceil(loc_cellcenter_i)), int(floor(loc_cellcenter_j))]
        grid_val_11 = self.arr[int(ceil(loc_cellcenter_i)), int(ceil(loc_cellcenter_j))]

        delta_i = loc_cellcenter_i - floor(loc_cellcenter_i)
        delta_j = loc_cellcenter_j - floor(loc_cellcenter_j)

        grid_val_y0 = grid_val_00 + (grid_val_10 - grid_val_00) * delta_i
        grid_val_y1 = grid_val_01 + (grid_val_11 - grid_val_01) * delta_i

        grid_val_interp = grid_val_y0 + (grid_val_y1 - grid_val_y0) * delta_j

        return grid_val_interp

    def magnitude_field(self):
        """
        Calculates magnitude field.

        :return: the magnitude fild
        """

        if not arr_check(self.grid_data):
            raise RasterIOExceptions("input requires data array with three dimensions and 2-level third dimension")

        vx = self.grid_data[:, :, 0]
        vy = self.grid_data[:, :, 1]

        magnitude = np.sqrt(vx ** 2 + vy ** 2)

        magnitude_fld = GeoArray()
        magnitude_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        magnitude_fld.set_grid_data(magnitude)

        return magnitude_fld

    # calculates orientations
    def orientations(self):

        if not arr_check(self.grid_data):
            raise RasterIOExceptions("input requires data array with three dimensions and 2-level third dimension")

        vx, vy = self.grid_data[:, :, 0], self.grid_data[:, :, 1]

        orientations = np.degrees(np.arctan2(vx, vy))
        orientations = np.where(orientations < 0.0, orientations + 360.0, orientations)

        orientations_fld = GeoArray()
        orientations_fld.set_grid_domain(self.get_grid_domain().get_start_point(),
                                         self.get_grid_domain().get_end_point())
        orientations_fld.set_grid_data(orientations)

        return orientations_fld

    # calculates divergence_2D
    def divergence(self):

        if not arr_check(self.grid_data):
            raise RasterIOExceptions("input requires data array with three dimensions and 2-level third dimension")

        dvx_dx = np.gradient(self.grid_data[:, :, 0])[1]
        dvy_dy = -(np.gradient(self.grid_data[:, :, 1])[0])

        divergence_fld = GeoArray()
        divergence_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        divergence_fld.set_grid_data((dvx_dx + dvy_dy) / self.get_cellsize_horiz_mean())

        return divergence_fld

    # calculates curl module
    def curl_module(self):

        if not arr_check(self.grid_data):
            raise RasterIOExceptions("input requires data array with three dimensions and 2-level third dimension")

        dvy_dx = np.gradient(self.grid_data[:, :, 1])[1]
        dvx_dy = -(np.gradient(self.grid_data[:, :, 0])[0])

        curl_fld = GeoArray()
        curl_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        curl_fld.set_grid_data((dvy_dx - dvx_dy) / self.get_cellsize_horiz_mean())

        return curl_fld

    # calculates magnitude gradient along x axis
    def grad_xaxis(self):

        if not arr_check(self.grid_data):
            raise RasterIOExceptions("input requires data array with three dimensions and 2-level third dimension")

        vx, vy = self.grid_data[:, :, 0], self.grid_data[:, :, 1]

        dir_array = np.arctan2(vx, vy)

        vect_magn = np.sqrt(vx ** 2 + vy ** 2)
        dm_dy, dm_dx = np.gradient(vect_magn)

        vect_xgrad_fld = GeoArray()
        vect_xgrad_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        vect_xgrad_fld.set_grid_data(dm_dx / self.get_cellsize_horiz_mean())

        return vect_xgrad_fld

    # calculates magnitude gradient along y axis
    def grad_yaxis(self):

        if not arr_check(self.grid_data):
            raise RasterIOExceptions("input requires data array with three dimensions and 2-level third dimension")

        vx, vy = self.grid_data[:, :, 0], self.grid_data[:, :, 1]

        dir_array = np.arctan2(vx, vy)

        vect_magn = np.sqrt(vx ** 2 + vy ** 2)
        dm_dy, dm_dx = np.gradient(vect_magn)
        dm_dy = - dm_dy

        vect_ygrad_fld = GeoArray()
        vect_ygrad_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        vect_ygrad_fld.set_grid_data(dm_dy / self.get_cellsize_horiz_mean())

        return vect_ygrad_fld

    # calculates gradient along flow lines
    def grad_flowlines(self):

        if not arr_check(self.grid_data):
            raise RasterIOExceptions("input requires data array with three dimensions and 2-level third dimension")

        vx, vy = self.grid_data[:, :, 0], self.grid_data[:, :, 1]

        dir_array = np.arctan2(vx, vy)

        vect_magn = np.sqrt(vx ** 2 + vy ** 2)
        dm_dy, dm_dx = np.gradient(vect_magn)
        dm_dy = - dm_dy

        velocity_gradient = dm_dx * np.sin(dir_array) + dm_dy * np.cos(dir_array)

        vect_magn_grad_fld = GeoArray()
        vect_magn_grad_fld.set_grid_domain(self.get_grid_domain().get_start_point(),
                                           self.get_grid_domain().get_end_point())
        vect_magn_grad_fld.set_grid_data(velocity_gradient / self.get_cellsize_horiz_mean())

        return vect_magn_grad_fld

    # returns the velocity components interpolated using the bilinear interpolation
    def interpolate_level_bilinear(self, level, curr_Pt_gridcoord):

        try:
            v_00 = self.get_grid_data(level)[floor(curr_Pt_gridcoord.i - 0.5), floor(curr_Pt_gridcoord.j - 0.5)]
            v_10 = self.get_grid_data(level)[floor(curr_Pt_gridcoord.i - 0.5), ceil(curr_Pt_gridcoord.j - 0.5)]
            v_01 = self.get_grid_data(level)[ceil(curr_Pt_gridcoord.i - 0.5), floor(curr_Pt_gridcoord.j - 0.5)]
            v_11 = self.get_grid_data(level)[ceil(curr_Pt_gridcoord.i - 0.5), ceil(curr_Pt_gridcoord.j - 0.5)]
        except:
            raise GeoArrayIOException('Error in get_grid_data() function')

        delta_j_grid = curr_Pt_gridcoord.j - int(curr_Pt_gridcoord.j)

        assert delta_j_grid >= 0 and delta_j_grid < 1

        if delta_j_grid >= 0.5:
            delta_j_grid = delta_j_grid - 0.5
        else:
            delta_j_grid = 0.5 + delta_j_grid

        v_y0 = v_00 + (v_10 - v_00) * delta_j_grid
        v_y1 = v_01 + (v_11 - v_01) * delta_j_grid

        delta_i_grid = curr_Pt_gridcoord.i - int(curr_Pt_gridcoord.i)

        assert delta_i_grid >= 0 and delta_i_grid < 1

        if delta_i_grid >= 0.5:
            delta_i_grid = delta_i_grid - 0.5
        else:
            delta_i_grid = 0.5 + delta_i_grid

        interp_v = v_y0 + (v_y1 - v_y0) * delta_i_grid

        return interp_v

        # check if a point pathline can be evaluated

    def include_point_location(self, currPoint):

        curr_Pt_gridcoord = self.geog2gridcoord(currPoint)

        if curr_Pt_gridcoord.i < 0.5 or curr_Pt_gridcoord.j < 0.5 or \
                curr_Pt_gridcoord.j > (self.get_xlines_num() - 0.5) or curr_Pt_gridcoord.i > (
                self.get_ylines_num() - 0.5):
            return False

        if len(self.grid_data.shape) != 3 and \
                self.grid_data.shape[2] != 2:
            raise RasterIOExceptions("input requires data array with three dimensions and 2-level third dimension")

        for n in range(2):
            if np.isnan(self.get_grid_data(n)[floor(curr_Pt_gridcoord.i - 0.5), floor(curr_Pt_gridcoord.j - 0.5)]) or \
                    np.isnan(
                        self.get_grid_data(n)[ceil(curr_Pt_gridcoord.i - 0.5), floor(curr_Pt_gridcoord.j - 0.5)]) or \
                    np.isnan(
                        self.get_grid_data(n)[floor(curr_Pt_gridcoord.i - 0.5), ceil(curr_Pt_gridcoord.j - 0.5)]) or \
                    np.isnan(self.get_grid_data(n)[ceil(curr_Pt_gridcoord.i - 0.5), ceil(curr_Pt_gridcoord.j - 0.5)]):
                return False

        return True

    def point_velocity(self, currpoint):
        """
        return the velocity components of a point in a velocity field
        """

        if not self.include_point_location(currpoint):
            return None, None

        currpoint_gridcoord = self.geog2gridcoord(currpoint)
        currpoint_vx = self.interpolate_level_bilinear(0, currpoint_gridcoord)
        currpoint_vy = self.interpolate_level_bilinear(1, currpoint_gridcoord)

        return currpoint_vx, currpoint_vy

    def interpolate_RKF(self, delta_time, curr_Pt):
        """
        interpolate points according to RKF method
        """

        K1_vx, K1_vy = self.point_velocity(curr_Pt)
        if K1_vx is None or K1_vy is None:
            return None, None

        K2_Pt = Point(curr_Pt.x + (0.25) * delta_time * K1_vx, curr_Pt.y + (0.25) * delta_time * K1_vy)
        K2_vx, K2_vy = self.point_velocity(K2_Pt)
        if K2_vx is None or K2_vy is None:
            return None, None

        K3_Pt = Point(curr_Pt.x + (3.0 / 32.0) * delta_time * K1_vx + (9.0 / 32.0) * delta_time * K2_vx, \
                      curr_Pt.y + (3.0 / 32.0) * delta_time * K1_vy + (9.0 / 32.0) * delta_time * K2_vy)
        K3_vx, K3_vy = self.point_velocity(K3_Pt)
        if K3_vx is None or K3_vy is None:
            return None, None

        K4_Pt = Point(curr_Pt.x + (1932.0 / 2197.0) * delta_time * K1_vx - (7200.0 / 2197.0) * delta_time * K2_vx + (
                    7296.0 / 2197.0) * delta_time * K3_vx, \
                      curr_Pt.y + (1932.0 / 2197.0) * delta_time * K1_vy - (7200.0 / 2197.0) * delta_time * K2_vy + (
                                  7296.0 / 2197.0) * delta_time * K3_vy)
        K4_vx, K4_vy = self.point_velocity(K4_Pt)
        if K4_vx is None or K4_vy is None:
            return None, None

        K5_Pt = Point(curr_Pt.x + (439.0 / 216.0) * delta_time * K1_vx - (8.0) * delta_time * K2_vx + (
                    3680.0 / 513.0) * delta_time * K3_vx - (845.0 / 4104.0) * delta_time * K4_vx, \
                      curr_Pt.y + (439.0 / 216.0) * delta_time * K1_vy - (8.0) * delta_time * K2_vy + (
                                  3680.0 / 513.0) * delta_time * K3_vy - (845.0 / 4104.0) * delta_time * K4_vy)
        K5_vx, K5_vy = self.point_velocity(K5_Pt)
        if K5_vx is None or K5_vy is None:
            return None, None

        K6_Pt = Point(curr_Pt.x - (8.0 / 27.0) * delta_time * K1_vx + (2.0) * delta_time * K2_vx - (
                    3544.0 / 2565.0) * delta_time * K3_vx + (1859.0 / 4104.0) * delta_time * K4_vx - (
                                  11.0 / 40.0) * delta_time * K5_vx, \
                      curr_Pt.y - (8.0 / 27.0) * delta_time * K1_vy + (2.0) * delta_time * K2_vy - (
                                  3544.0 / 2565.0) * delta_time * K3_vy + (1859.0 / 4104.0) * delta_time * K4_vy - (
                                  11.0 / 40.0) * delta_time * K5_vy)
        K6_vx, K6_vy = self.point_velocity(K6_Pt)
        if K6_vx is None or K6_vy is None:
            return None, None

        rkf_4o_x = curr_Pt.x + delta_time * (
                    (25.0 / 216.0) * K1_vx + (1408.0 / 2565.0) * K3_vx + (2197.0 / 4104.0) * K4_vx - (
                        1.0 / 5.0) * K5_vx)
        rkf_4o_y = curr_Pt.y + delta_time * (
                    (25.0 / 216.0) * K1_vy + (1408.0 / 2565.0) * K3_vy + (2197.0 / 4104.0) * K4_vy - (
                        1.0 / 5.0) * K5_vy)
        temp_Pt = Point(rkf_4o_x, rkf_4o_y)

        interp_Pt_x = curr_Pt.x + delta_time * (
                    (16.0 / 135.0) * K1_vx + (6656.0 / 12825.0) * K3_vx + (28561.0 / 56430.0) * K4_vx - (
                        9.0 / 50.0) * K5_vx + (2.0 / 55.0) * K6_vx)
        interp_Pt_y = curr_Pt.y + delta_time * (
                    (16.0 / 135.0) * K1_vy + (6656.0 / 12825.0) * K3_vy + (28561.0 / 56430.0) * K4_vy - (
                        9.0 / 50.0) * K5_vy + (2.0 / 55.0) * K6_vy)
        interp_Pt = Point(interp_Pt_x, interp_Pt_y)

        interp_PT_error_estimate = interp_Pt.distance(temp_Pt)

        return interp_Pt, interp_PT_error_estimate

    # writes ESRI ascii grid
    def write_esrigrid(self, outgrid_fn, esri_nullvalue=-99999, level=0):

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
        outputgrid.write("NCOLS %d\n" % self.get_xlines_num())
        outputgrid.write("NROWS %d\n" % self.get_ylines_num())
        outputgrid.write("XLLCORNER %.8f\n" % self.grid_domain.get_start_point().x)
        outputgrid.write("YLLCORNER %.8f\n" % self.grid_domain.get_start_point().y)
        outputgrid.write("CELLSIZE %.8f\n" % self.get_cellsize_horiz_mean())
        outputgrid.write("NODATA_VALUE %d\n" % esri_nullvalue)

        esrigrid_outvalues = np.where(np.isnan(self.grid_data), esri_nullvalue, self.grid_data)

        # output of results
        if len(self.grid_data.shape) == 3:
            for i in range(0, self.get_ylines_num()):
                for j in range(0, self.get_xlines_num()):
                    outputgrid.write("%.8f " % (esrigrid_outvalues[i, j, level]))
                outputgrid.write("\n")
        elif len(self.grid_data.shape) == 2:
            for i in range(0, self.get_ylines_num()):
                for j in range(0, self.get_xlines_num()):
                    outputgrid.write("%.8f " % (esrigrid_outvalues[i, j]))
                outputgrid.write("\n")

        outputgrid.close()

        return True


if __name__ == "__main__":

    import doctest
    doctest.testmod()
