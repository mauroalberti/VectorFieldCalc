
from math import ceil, floor
import copy

import numpy as np

from ..gsf.geometry import MIN_SEPARATION_THRESHOLD, Point


class ArrCoord(object):
    """
    2D Array coordinates.
    Manages coordinates in the raster (array) space.
    """

    def __init__(self, ival=0.0, jval=0.0):
        """
        :param  ival:  the i (-y) array coordinate of the point.
        :type  ival:  number or string convertible to float.
        :param  jval:  the j (x) array coordinate of the point.
        :type  jval:  number or string convertible to float.

        :return:  self.
        """
        self._i = float(ival)
        self._j = float(jval)

    def g_i(self):
        """
        Get i (row) coordinate value.

        :return:  the i (-y) array coordinate of the point - float.
        """
        return self._i

    def s_i(self, ival):
        """
        Set i (row) coordinate value.

        :param  ival:  the i (-y) array coordinate of the point.
        :type  ival:  number or string convertible to float.

        :return:  self.
        """
        self._i = float(ival)

    # set property for i
    i = property(g_i, s_i)

    def g_j(self):
        """
        Get j (column) coordinate value.

        :return:  the j (x) array coordinate of the point - float.
        """
        
        return self._j

    def s_j(self, jval):
        """
        Set j (column) coordinate value.

        :param  jval:  the j (x) array coordinate of the point.
        :type  jval:  number or string convertible to float.

        :return:  self.
        """
        
        self._j = jval

    j = property(g_j, s_j)

    def grid2geogcoord(self, currGeoGrid):
        """
        
        :param currGeoGrid: 
        :return: 
        """
        
        currPt_geogr_y = currGeoGrid.domain.trcorner.p_y - self.i * currGeoGrid.cellsize_y
        currPt_geogr_x = currGeoGrid.domain.llcorner.p_x + self.j * currGeoGrid.cellsize_x

        return currPt_geogr_x, currPt_geogr_y


class RectangularDomain(object):
    """
    Rectangular spatial domain class.

    """

    def __init__(self, pt_llc=None, pt_trc=None):
        """
        Class constructor.

        :param  pt_llc:  lower-left corner of the domain.
        :type  pt_llc:  Point.
        :param  pt_trc:  top-right corner of the domain.
        :type  pt_trc:  Point.

        :return:  RectangularDomain instance.
        """
        
        self._llcorner = pt_llc
        self._trcorner = pt_trc

    @property
    def llcorner(self):
        """
        Get lower-left corner of the spatial domain.

        :return:  lower-left corner of the spatial domain - Point.
        """
        
        return self._llcorner

    @property
    def trcorner(self):
        """
        Get top-right corner of the spatial domain.

        :return:  top-right corner of the spatial domain - Point.
        """
        
        return self._trcorner

    @property
    def xrange(self):
        """
        Get x range of spatial domain.

        :return:  x range - float.
        """
        
        return self.trcorner.p_x - self.llcorner.p_x

    @property
    def yrange(self):
        """
        Get y range of spatial domain.

        :return:  y range - float.
        """
        
        return self.trcorner.p_y - self.llcorner.p_y

    @property
    def zrange(self):
        """
        Get z range of spatial domain.

        :return:  z range - float.
        """
        
        return self.trcorner.p_z - self.llcorner.p_z

    @property
    def horiz_area(self):
        """
        Get horizontal area of spatial domain.

        :return:  area - float.
        """

        return self.xrange * self.yrange


def arr_check(curr_array):
    """
    check array input.

    :param curr_array:
    :return:
    """
    if len(curr_array.shape) != 3 and \
            curr_array.shape[2] != 2:
        return False
    return True


class Grid(object):
    """
    Grid class.
    Stores and manages the most of data and processing.

    """

    def __init__(self, source_filename=None, grid_params=None, grid_data=None):
        """
        Grid class constructor.

        :param  source_filename:  name of file from which data and geo-parameters derive.
        :type  source_filename:  string.
        :param  grid_params:  the geo-parameters of the grid.
        :type  grid_params:  class GDALParameters.
        :param  grid_data:  the array storing the data.
        :type  grid_data:  2D np.array.

        :return:  self.
        """
        self._sourcename = source_filename

        if grid_params is not None:
            pt_llc = grid_params.llcorner
            pt_trc = grid_params.trcorner
        else:
            pt_llc = None
            pt_trc = None

        self._grid_domain = RectangularDomain(pt_llc, pt_trc)

        if grid_data is not None:
            self._grid_data = grid_data.copy()
        else:
            self._grid_data = None

    def s_domain(self, domain):
        """
        Set spatial domain.

        :param  domain:  Spatial domain to be attributed to the current Grid instance.
        :type  domain:  class RectangularDomain.

        :return: self
        """

        del self._grid_domain
        self._grid_domain = copy.deepcopy(domain)

    def g_domain(self):
        """
        Get spatial domain.

        :return: the spatial domain of the current Grid instance - class RectangularDomain.
        """

        return self._grid_domain

    def d_domain(self):
        """
        Delete current spatial domain of the Grid instance.

        :return: self
        """

        del self._grid_domain

    # set property for spatial domain
    domain = property(g_domain, s_domain, d_domain)

    def s_grid_data(self, data_array):
        """
        Set grid data array.

        :param data_array: numpy.array of data values.
        :param type: 2D numpy.array.

        :return: self.
        """

        if self._grid_data is not None:
            del self._grid_data

        self._grid_data = data_array.copy()

    def g_grid_data(self):
        """
        Get grid data array.

        :return: 2D numpy.array.
        """

        return self._grid_data

    def d_grid_data(self):
        """
        Delete grid data array.

        :return: self.
        """

        del self._grid_data

    data = property(g_grid_data, s_grid_data, d_grid_data)

    def grid_extent(self):
        """
        Return the xmin, xmax and ymin, ymax values as a dictionary
        """

        return dict(xmin=self.domain.llcorner.p_x,
                    xmax=self.domain.trcorner.p_x,
                    ymin=self.domain.llcorner.p_y,
                    ymax=self.domain.trcorner.p_y)

    @property
    def xmin(self):

        return self.grid_extent()['xmin']

    @property
    def xmax(self):

        return self.grid_extent()['xmax']

    @property
    def ymin(self):

        return self.grid_extent()['ymin']

    @property
    def ymax(self):

        return self.grid_extent()['ymax']

    @property
    def row_num(self):
        """
        Get row number of the grid domain.

        :return: number of rows of data array - int.
        """

        return np.shape(self.data)[0]

    @property
    def col_num(self):
        """
        Get column number of the grid domain.

        :return: number of columns of data array - int.
        """

        return np.shape(self.data)[1]

    @property
    def cellsize_x(self):
        """
        Get the cell size of the grid in the x direction.

        :return: cell size in the x (j) direction - float.
        """

        return self.domain.xrange / float(self.col_num)

    @property
    def cellsize_y(self):
        """
        Get the cell size of the grid in the y direction.

        :return: cell size in the y (-i) direction - float.
        """

        return self.domain.yrange / float(self.row_num)

    @property
    def cellsize_h(self):
        """
        Get the mean horizontal cell size.

        :return: mean horizontal cell size - float.
        """

        return (self.cellsize_x + self.cellsize_y) / 2.0

    def geog2array_coord(self, curr_Pt):
        """
        Converts from geographic to raster (array) coordinates.

        :param curr_Pt: point whose geographical coordinates will be converted to raster (array) ones.
        :type curr_Pt: Point.

        :return: point coordinates in raster (array) frame - class ArrCoord.
        """
        currArrCoord_grid_i = (self.domain.trcorner.p_y - curr_Pt.p_y) / self.cellsize_y
        currArrCoord_grid_j = (curr_Pt.p_x - self.domain.llcorner.p_x) / self.cellsize_x

        return ArrCoord(currArrCoord_grid_i, currArrCoord_grid_j)

    def x(self):
        """
        Creates an array storing the geographical coordinates of the cell centers along the x axis.
        Direction is from left to right.

        :return: numpy.array, shape: 1 x col_num.
        """

        x_values = self.domain.llcorner.p_x + self.cellsize_x * (0.5 + np.arange(self.col_num))

        return x_values[np.newaxis, :]

    def y(self):
        """
        Creates an array storing the geographical coordinates of the cell centers along the y axis.
        Direction is from top to bottom.

        :return: numpy.array, shape: row_num x 1.
        """

        y_values = self.domain.trcorner.p_y - self.cellsize_y * (0.5 + np.arange(self.row_num))

        return y_values[:, np.newaxis]

    def grad_forward_y(self):
        """
        Return an array representing the forward gradient in the y direction (top-wards), with values scaled by cell size.

        :return: numpy.array, same shape as current Grid instance
        """

        gf = np.zeros(np.shape(self.data)) * np.NaN
        gf[1:, :] = self.data[:-1, :] - self.data[1:, :]

        return gf / float(self.cellsize_y)

    def grad_forward_x(self):
        """
        Return an array representing the forward gradient in the x direction (right-wards), with values scaled by cell size.

        :return: numpy.array, same shape as current Grid instance
        """

        gf = np.zeros(np.shape(self.data), ) * np.NaN
        gf[:, :-1] = self.data[:, 1:] - self.data[:, :-1]

        return gf / float(self.cellsize_x)

    def interpolate_bilinear(self, curr_Pt_array_coord):
        """
        Interpolate the z value at a point, given its array coordinates.
        Interpolation method: bilinear.

        :param curr_Pt_array_coord: array coordinates of the point for which the interpolation will be made.
        :type curr_Pt_array_coord: class ArrCoord.

        :return: interpolated z value - float.
        """

        currPt_cellcenter_i = curr_Pt_array_coord.i - 0.5
        currPt_cellcenter_j = curr_Pt_array_coord.j - 0.5

        assert currPt_cellcenter_i > 0, currPt_cellcenter_j > 0

        grid_val_00 = self.data[int(floor(currPt_cellcenter_i)), int(floor(currPt_cellcenter_j))]
        grid_val_01 = self.data[int(floor(currPt_cellcenter_i)), int(ceil(currPt_cellcenter_j))]
        grid_val_10 = self.data[int(ceil(currPt_cellcenter_i)), int(floor(currPt_cellcenter_j))]
        grid_val_11 = self.data[int(ceil(currPt_cellcenter_i)), int(ceil(currPt_cellcenter_j))]

        delta_i = currPt_cellcenter_i - floor(currPt_cellcenter_i)
        delta_j = currPt_cellcenter_j - floor(currPt_cellcenter_j)

        grid_val_y0 = grid_val_00 + (grid_val_10 - grid_val_00) * delta_i
        grid_val_y1 = grid_val_01 + (grid_val_11 - grid_val_01) * delta_i

        grid_val_interp = grid_val_y0 + (grid_val_y1 - grid_val_y0) * delta_j

        return grid_val_interp

    def intersection_with_surface(self, surf_type, srcPt, srcPlaneAttitude):
        """
        Calculates the intersections (as points) between DEM (the self object) and an analytical surface.
        Currently it works only with planes.

        :param surf_type: type of considered surface (i.e., plane, the only case implemented at present).
        :type surf_type: String.
        :param srcPt: point, expressed in geographical coordinates, that the plane must contain.
        :type srcPt: Point.
        :param srcPlaneAttitude: orientation of the surface (currently only planes).
        :type srcPlaneAttitude: class GPlane.

        :return: tuple of four arrays
        """

        if surf_type == 'plane':

            # closures to compute the geographic coordinates (in x- and y-) of a cell center
            # the grid coordinates of the cell center are expressed by i and j
            grid_coord_to_geogr_coord_x_closure = lambda j: self.domain.llcorner.p_x + self.cellsize_x * (0.5 + j)
            grid_coord_to_geogr_coord_y_closure = lambda i: self.domain.trcorner.p_y - self.cellsize_y * (0.5 + i)

            # arrays storing the geographical coordinates of the cell centers along the x- and y- axes
            cell_center_x_array = self.x()
            cell_center_y_array = self.y()

            ycoords_x, xcoords_y = np.broadcast_arrays(cell_center_x_array, cell_center_y_array)

            #### x-axis direction intersections

            # 2D array of DEM segment parameters
            x_dem_m = self.grad_forward_x()
            x_dem_q = self.data - cell_center_x_array * x_dem_m

            # closure for the planar surface that, given (x,y), will be used to derive z
            plane_z_closure = srcPlaneAttitude.plane_from_geo(srcPt)

            # 2D array of plane segment parameters
            x_plane_m = srcPlaneAttitude.plane_x_coeff()
            x_plane_q = np.array_from_function(self.row_num(), 1, lambda j: 0, grid_coord_to_geogr_coord_y_closure,
                                               plane_z_closure)

            # 2D array that defines denominator for intersections between local segments
            x_inters_denomin = np.where(x_dem_m != x_plane_m, x_dem_m - x_plane_m, np.NaN)

            coincident_x = np.where(x_dem_q != x_plane_q, np.NaN, ycoords_x)

            xcoords_x = np.where(x_dem_m != x_plane_m, (x_plane_q - x_dem_q) / x_inters_denomin, coincident_x)
            xcoords_x = np.where(xcoords_x < ycoords_x, np.NaN, xcoords_x)
            xcoords_x = np.where(xcoords_x >= ycoords_x + self.cellsize_x, np.NaN, xcoords_x)

            #### y-axis direction intersections

            # 2D array of DEM segment parameters
            y_dem_m = self.grad_forward_y()
            y_dem_q = self.data - cell_center_y_array * y_dem_m

            # 2D array of plane segment parameters
            y_plane_m = srcPlaneAttitude.plane_y_coeff()
            y_plane_q = np.array_from_function(1, self.col_num, grid_coord_to_geogr_coord_x_closure, lambda i: 0,
                                               plane_z_closure)

            # 2D array that defines denominator for intersections between local segments
            y_inters_denomin = np.where(y_dem_m != y_plane_m, y_dem_m - y_plane_m, np.NaN)
            coincident_y = np.where(y_dem_q != y_plane_q, np.NaN, xcoords_y)

            ycoords_y = np.where(y_dem_m != y_plane_m, (y_plane_q - y_dem_q) / y_inters_denomin, coincident_y)

            # filter out cases where intersection is outside cell range
            ycoords_y = np.where(ycoords_y < xcoords_y, np.NaN, ycoords_y)
            ycoords_y = np.where(ycoords_y >= xcoords_y + self.cellsize_y, np.NaN, ycoords_y)

            for i in xrange(xcoords_x.shape[0]):
                for j in xrange(xcoords_x.shape[1]):
                    if abs(xcoords_x[i, j] - ycoords_x[i, j]) < MIN_SEPARATION_THRESHOLD and abs(
                                    ycoords_y[i, j] - xcoords_y[i, j]) < MIN_SEPARATION_THRESHOLD:
                        ycoords_y[i, j] = np.NaN

            return xcoords_x, xcoords_y, ycoords_x, ycoords_y


    # calculates magnitude
    def magnitude(self):

        if not arr_check(self.grid_data):
            raise FunInp_Err, 'input requires data array with three dimensions and 2-level third dimension'

        vx = self.grid_data[:, :, 0]
        vy = self.grid_data[:, :, 1]

        magnitude = np.sqrt(vx ** 2 + vy ** 2)

        magnitude_fld = Grid()
        magnitude_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        magnitude_fld.set_grid_data(magnitude)

        return magnitude_fld

    # calculates orientations
    def orientations(self):

        if not arr_check(self.grid_data):
            raise FunInp_Err, 'input requires data array with three dimensions and 2-level third dimension'

        vx, vy = self.grid_data[:, :, 0], self.grid_data[:, :, 1]

        orientations = np.degrees(np.arctan2(vx, vy))
        orientations = np.where(orientations < 0.0, orientations + 360.0, orientations)

        orientations_fld = Grid()
        orientations_fld.set_grid_domain(self.get_grid_domain().get_start_point(),
                                         self.get_grid_domain().get_end_point())
        orientations_fld.set_grid_data(orientations)

        return orientations_fld

    # calculates divergence
    def divergence(self):

        if not arr_check(self.grid_data):
            raise FunInp_Err, 'input requires data array with three dimensions and 2-level third dimension'

        dvx_dx = np.gradient(self.grid_data[:, :, 0])[1]
        dvy_dy = -(np.gradient(self.grid_data[:, :, 1])[0])

        divergence_fld = Grid()
        divergence_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        divergence_fld.set_grid_data((dvx_dx + dvy_dy) / self.get_cellsize_horiz_mean())

        return divergence_fld

    # calculates curl module
    def curl_module(self):

        if not arr_check(self.grid_data):
            raise FunInp_Err, 'input requires data array with three dimensions and 2-level third dimension'

        dvy_dx = np.gradient(self.grid_data[:, :, 1])[1]
        dvx_dy = -(np.gradient(self.grid_data[:, :, 0])[0])

        curl_fld = Grid()
        curl_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        curl_fld.set_grid_data((dvy_dx - dvx_dy) / self.get_cellsize_horiz_mean())

        return curl_fld

    # calculates magnitude gradient along x axis
    def grad_xaxis(self):

        if not arr_check(self.grid_data):
            raise FunInp_Err, 'input requires data array with three dimensions and 2-level third dimension'

        vx, vy = self.grid_data[:, :, 0], self.grid_data[:, :, 1]

        dir_array = np.arctan2(vx, vy)

        vect_magn = np.sqrt(vx ** 2 + vy ** 2)
        dm_dy, dm_dx = np.gradient(vect_magn)

        vect_xgrad_fld = Grid()
        vect_xgrad_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        vect_xgrad_fld.set_grid_data(dm_dx / self.get_cellsize_horiz_mean())

        return vect_xgrad_fld

    # calculates magnitude gradient along y axis
    def grad_yaxis(self):

        if not arr_check(self.grid_data):
            raise FunInp_Err, 'input requires data array with three dimensions and 2-level third dimension'

        vx, vy = self.grid_data[:, :, 0], self.grid_data[:, :, 1]

        dir_array = np.arctan2(vx, vy)

        vect_magn = np.sqrt(vx ** 2 + vy ** 2)
        dm_dy, dm_dx = np.gradient(vect_magn)
        dm_dy = - dm_dy

        vect_ygrad_fld = Grid()
        vect_ygrad_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        vect_ygrad_fld.set_grid_data(dm_dy / self.get_cellsize_horiz_mean())

        return vect_ygrad_fld

    # calculates gradient along flow lines
    def grad_flowlines(self):

        if not arr_check(self.grid_data):
            raise FunInp_Err, 'input requires data array with three dimensions and 2-level third dimension'

        vx, vy = self.grid_data[:, :, 0], self.grid_data[:, :, 1]

        dir_array = np.arctan2(vx, vy)

        vect_magn = np.sqrt(vx ** 2 + vy ** 2)
        dm_dy, dm_dx = np.gradient(vect_magn)
        dm_dy = - dm_dy

        velocity_gradient = dm_dx * np.sin(dir_array) + dm_dy * np.cos(dir_array)

        vect_magn_grad_fld = Grid()
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
            raise Raster_Parameters_Errors, 'Error in get_grid_data() function'

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
            raise FunInp_Err, 'input requires data array with three dimensions and 2-level third dimension'

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
            raise Output_Errors, "Output grid '%s' already exists" % outgrid_fn

        try:
            outputgrid = open(outgrid_fn, 'w')  # create the output ascii file
        except:
            raise Output_Errors, "Unable to create output grid '%s'" % outgrid_fn

        if outputgrid is None:
            raise Output_Errors, "Unable to create output grid '%s'" % outgrid_fn

        # writes header of grid ascii file
        outputgrid.write('NCOLS %d\n' % self.get_xlines_num())
        outputgrid.write('NROWS %d\n' % self.get_ylines_num())
        outputgrid.write('XLLCORNER %.8f\n' % self.grid_domain.get_start_point().x)
        outputgrid.write('YLLCORNER %.8f\n' % self.grid_domain.get_start_point().y)
        outputgrid.write('CELLSIZE %.8f\n' % self.get_cellsize_horiz_mean())
        outputgrid.write('NODATA_VALUE %d\n' % esri_nullvalue)

        esrigrid_outvalues = np.where(np.isnan(self.grid_data), esri_nullvalue, self.grid_data)

        # output of results
        if len(self.grid_data.shape) == 3:
            for i in range(0, self.get_ylines_num()):
                for j in range(0, self.get_xlines_num()):
                    outputgrid.write('%.8f ' % (esrigrid_outvalues[i, j, level]))
                outputgrid.write('\n')
        elif len(self.grid_data.shape) == 2:
            for i in range(0, self.get_ylines_num()):
                for j in range(0, self.get_xlines_num()):
                    outputgrid.write('%.8f ' % (esrigrid_outvalues[i, j]))
                outputgrid.write('\n')

        outputgrid.close()

        return True
