# -*- coding: utf-8 -*-


import copy

from ...mathematics.defaults import *
from ...defaults.typing import *

from .exceptions import *

from ..generics.general import RectangularDomain


class ArrCoord(object):
    """
    2D Array coordinates.
    Manages coordinates in the rasters (array) space.
    """

    def __init__(self, ival: float=0.0, jval:float=0.0) -> None:
        """
        :param  ival:  the i (-y) array coordinate of the point.
        :type  ival:  number or string convertible to float.
        :param  jval:  the j (x) array coordinate of the point.
        :type  jval:  number or string convertible to float.

        :return:  self.
        """
        self._i = float(ival)
        self._j = float(jval)

    def g_i(self) -> float:
        """
        Get i (row) coordinate value.

        :return:  the i (-y) array coordinate of the point - float.
        """
        return self._i

    def s_i(self, ival) -> None:
        """
        Set i (row) coordinate value.

        :param  ival:  the i (-y) array coordinate of the point.
        :type  ival:  number or string convertible to float.

        :return:  self.
        """
        self._i = float(ival)

    # set property for i
    i = property(g_i, s_i)

    def g_j(self) -> float:
        """
        Get j (column) coordinate value.

        :return:  the j (x) array coordinate of the point - float.
        """
        
        return self._j

    def s_j(self, jval) -> None:
        """
        Set j (column) coordinate value.

        :param  jval:  the j (x) array coordinate of the point.
        :type  jval:  number or string convertible to float.

        :return:  self.
        """
        
        self._j = jval

    j = property(g_j, s_j)

    def grid2geogcoord(self, currGeoGrid : 'Grid') -> Tuple[float, float]:
        """
        
        :param currGeoGrid: a Grid instance.
        :return: the x and y locations.
        :rtype: tuple with two floats.
        """
        
        currPt_geogr_y = currGeoGrid.domain.trcorner.y - self.i * currGeoGrid.cellsize_y
        currPt_geogr_x = currGeoGrid.domain.llcorner.x + self.j * currGeoGrid.cellsize_x

        return currPt_geogr_x, currPt_geogr_y


def arr_check(curr_array):
    """
    Check array input.

    :param curr_array:
    :return:
    """
    if len(curr_array.shape) != 3 and \
            curr_array.shape[2] != 2:
        return False
    return True


class GridParams(object):
    """
    Class storing the grid parameters.
    """

    def __init__(
            self,
            inLlxy: Tuple[Number, Number],
            inUrxy: Optional[Tuple[Number, Number]]=None,
            inCellSizeX: Optional[Number]=None,
            inCellSizeY: Optional[Number]=None,
            inNumRows: Optional[int]=None,
            inNumCols: Optional[int]=None,
            inRotX: Optional[Number]=None,
            inRotY: Optional[Number]=None,
            inCrs: Optional[str]=None):
        """
        Creates the grid parameters instance.
        """

        self.llxy = inLlxy
        self.urxy = inUrxy
        self.cellsize_x = inCellSizeX
        self.cellsize_y = inCellSizeY
        self.num_rows = inNumRows
        self.num_cols = inNumCols
        self.rot_x = inRotX
        self.rot_y = inRotY
        self.crs = inCrs


class Grid(object):
    """
    Grid class.
    Stores and manages the most of data and processing.

    """

    def __init__(self, inGridParams: GridParams, inLevels: List['array']) -> None:
        """
        Grid class constructor.

        :param  inGridParams:  the geo-parameters of the grid.
        :type  inGridParams:  GridParams.
        :param  inLevels:  the list sotring the array levels.
        :type  inLevels:  list of 2D np.array.

        :return:  None.

        Examples:
        """

        if inGridParams is not None:
            pt_llc = inGridParams.llcorner
            pt_trc = inGridParams.trcorner
        else:
            pt_llc = None
            pt_trc = None

        self._grid_domain = RectangularDomain(pt_llc, pt_trc)

        if inLevels is not None:
            self._grid_levels = inLevels
        else:
            self._grid_levels = None

    def s_domain(self, domain: RectangularDomain) -> None:
        """
        Set spatial domain.

        :param  domain:  Spatial domain to be attributed to the current Grid instance.
        :type  domain: RectangularDomain.

        :return: None

        Examples:
        """

        del self._grid_domain
        self._grid_domain = copy.deepcopy(domain)

    def g_domain(self) -> RectangularDomain:
        """
        Get spatial domain.

        :return: the spatial domain of the current Grid instance
        ;:type: RectangularDomain.

        Examples:
        """

        return self._grid_domain

    def d_domain(self) -> None:
        """
        Delete current spatial domain of the Grid instance.

        :return: None

        Examples:
        """

        del self._grid_domain

    # set property for spatial domain
    domain = property(g_domain, s_domain, d_domain)

    def s_grid_data(self, data_array: 'array') -> None:
        """
        Set grid data array.

        :param data_array: numpy.array of data values.
        :param type: 2D numpy.array.

        :return: None

        Examples:
        """

        if self._grid_levels is not None:
            del self._grid_levels

        self._grid_levels = data_array.copy()

    def g_grid_data(self) -> 'array':
        """
        Get grid data array.

        :return: 2D array.
        :type: numpy.array

        Examples:
        """

        return self._grid_levels

    def d_grid_data(self) -> None:
        """
        Delete grid data array.

        :return: None.

        Examples:
        """

        del self._grid_levels

    data = property(g_grid_data, s_grid_data, d_grid_data)

    def grid_extent(self) -> Dict[str, float]:
        """
        Return the xmin, xmax and ymin, ymax values as a dictionary.

        :return: dictionary of values.
        :rtype: dictionary.
        """

        return dict(xmin=self.domain.llcorner.x,
                    xmax=self.domain.trcorner.x,
                    ymin=self.domain.llcorner.y,
                    ymax=self.domain.trcorner.y)

    @property
    def xmin(self) -> float:
        """
        Returns the x minimum value.

        :return:

        Examples:
        """

        return self.grid_extent()["xmin"]

    @property
    def xmax(self) -> float:
        """
        Returns the x maximum value.

        :return:

        Examples:
        """

        return self.grid_extent()["xmax"]

    @property
    def ymin(self) -> float:
        """
        Returns the y minimum value.

        :return:

        Examples:
        """

        return self.grid_extent()["ymin"]

    @property
    def ymax(self) -> float:
        """
        Returns the y maximum value.

        :return:

        Examples:
        """

        return self.grid_extent()["ymax"]

    @property
    def row_num(self) -> int:
        """
        Get row number of the grid domain.

        :return: number of rows of data array.
        :rtype: int.

        Examples:
        """

        return np.shape(self.data)[0]

    @property
    def col_num(self) -> int:
        """
        Get column number of the grid domain.

        :return: number of columns of data array
        :rtype: int.

        Examples:
        """

        return np.shape(self.data)[1]

    @property
    def cellsize_x(self) -> float:
        """
        Get the cell size of the grid in the x direction.

        :return: cell size in the x (j) direction
        :rtype: float.

        Examples:
        """

        return self.domain.xrange / float(self.col_num)

    @property
    def cellsize_y(self) -> float:
        """
        Get the cell size of the grid in the y direction.

        :return: cell size in the y (-i) direction
        :rtype: float.

        Examples:
        """

        return self.domain.yrange / float(self.row_num)

    @property
    def cellsize_h(self) -> float:
        """
        Get the mean horizontal cell size.

        :return: mean horizontal cell size
        :rtype: float.

        Examples:
        """

        return (self.cellsize_x + self.cellsize_y) / 2.0

    def geog2array_coord(self, x: float, y: float) -> ArrCoord:
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

        currArrCoord_grid_i = (self.domain.trcorner.y - y) / self.cellsize_y
        currArrCoord_grid_j = (x - self.domain.llcorner.x) / self.cellsize_x

        return ArrCoord(currArrCoord_grid_i, currArrCoord_grid_j)

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
          >>> Grid(grid_data=array([[1, 2], [3, 4]])).levels_num
          1
          >>> Grid(grid_data=np.ones((4, 3, 2)))
          2
        """

        print(type(self.data))
        return self.data.ndim

    def grad_forward_y(self) -> 'array':
        """
        Return an array representing the forward gradient in the y direction (top-wards), with values scaled by cell size.

        :return: array with same shape as current Grid instance
        :rtype: numpy.array

        Examples:
        """

        gf = np.zeros(np.shape(self.data)) * np.NaN
        gf[1:, :] = self.data[:-1, :] - self.data[1:, :]

        return gf / float(self.cellsize_y)

    def grad_forward_x(self) -> 'array':
        """
        Return an array representing the forward gradient in the x direction (right-wards), with values scaled by cell size.

        :return: array with same shape as current Grid instance
        :rtype: numpy.array

        Examples:
        """

        gf = np.zeros(np.shape(self.data), ) * np.NaN
        gf[:, :-1] = self.data[:, 1:] - self.data[:, :-1]

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

        grid_val_00 = self.data[int(floor(loc_cellcenter_i)), int(floor(loc_cellcenter_j))]
        grid_val_01 = self.data[int(floor(loc_cellcenter_i)), int(ceil(loc_cellcenter_j))]
        grid_val_10 = self.data[int(ceil(loc_cellcenter_i)), int(floor(loc_cellcenter_j))]
        grid_val_11 = self.data[int(ceil(loc_cellcenter_i)), int(ceil(loc_cellcenter_j))]

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

        magnitude_fld = Grid()
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

        orientations_fld = Grid()
        orientations_fld.set_grid_domain(self.get_grid_domain().get_start_point(),
                                         self.get_grid_domain().get_end_point())
        orientations_fld.set_grid_data(orientations)

        return orientations_fld

    # calculates divergence
    def divergence(self):

        if not arr_check(self.grid_data):
            raise RasterIOExceptions("input requires data array with three dimensions and 2-level third dimension")

        dvx_dx = np.gradient(self.grid_data[:, :, 0])[1]
        dvy_dy = -(np.gradient(self.grid_data[:, :, 1])[0])

        divergence_fld = Grid()
        divergence_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        divergence_fld.set_grid_data((dvx_dx + dvy_dy) / self.get_cellsize_horiz_mean())

        return divergence_fld

    # calculates curl module
    def curl_module(self):

        if not arr_check(self.grid_data):
            raise RasterIOExceptions("input requires data array with three dimensions and 2-level third dimension")

        dvy_dx = np.gradient(self.grid_data[:, :, 1])[1]
        dvx_dy = -(np.gradient(self.grid_data[:, :, 0])[0])

        curl_fld = Grid()
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

        vect_xgrad_fld = Grid()
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

        vect_ygrad_fld = Grid()
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
            raise RasterParametersExceptions('Error in get_grid_data() function')

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
