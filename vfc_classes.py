

from builtins import str
from builtins import range
from builtins import object
import os, sys, string
from math import *
import numpy as np
import gdal
from gdalconst import *


from .pygsf.libs_utils.gdal.gdal import *



# check if number:

def is_number(s):

    try:
        float(s)
        return True
    except:
        return False


# 2D Array coordinates 

class ArrCoord(object):
    
    # class constructor
    def __init__(self, ival, jval):
        self.i = ival
        self.j = jval

    # get i value
    def get_i(self):
        return self.i
   
    # get j value
    def get_j(self):
        return self.j
        
        
# 3D point class

class Point(object):
    
    # class constructor
    def __init__(self, x, y, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    # calculate euclidean distance between two points
    def distance(self, pt):
        return sqrt((self.x - pt.x)**2 + (self.y - pt.y)**2 + (self.z - pt.z)**2)
    
    # create a new point shifted by given amount
    def movedby(self, sx, sy, sz=0.0):
        return Point(self.x + sx, self.y + sy, self.z + sz)


# GDAL raster parameters
class GDALParameters(object):

    # class constructor
    def __init__(self): 
        self.driver_shortname = None
        self.datatype = None
        self.nodatavalue = None
        self.projection = ''
        self.topleftX = None
        self.topleftY = None
        self.pixsizeEW = None
        self.pixsizeNS = None
        self.rows = None
        self.cols = None
        self.rotationA = 0.0
        self.rotationB = 0.0
    
    # set driver format
    def set_driverShortName(self, raster_driver):
        self.driver_shortname = raster_driver
   
    # get driver format
    def get_driverShortName(self):
        return self.driver_shortname

    # set data type
    def set_dataType(self, data_type):
        self.datatype = data_type
   
    # get data type
    def get_dataType(self):
        return self.datatype

    # set nodata value
    def set_noDataValue(self, nodataval):
        self.nodatavalue = nodataval
   
    # get nodata value
    def get_noDataValue(self):
        return self.nodatavalue

    # set projection
    def set_projection(self, proj):
        self.projection = proj
   
    # get projection
    def get_projection(self):
        return self.projection

    # set topleftX value
    def set_topLeftX(self, topleftX):
        self.topleftX = topleftX
   
    # get topleftX value
    def get_topLeftX(self):
        return self.topleftX

    # set topleftY value
    def set_topLeftY(self, topleftY):
        self.topleftY = topleftY
   
    # get topleftY value
    def get_topLeftY(self):
        return self.topleftY

    # set pixsizeEW value
    def set_pixSizeEW(self, pixsizeEW):
        self.pixsizeEW = pixsizeEW
   
    # get pixsizeEW value
    def get_pixSizeEW(self):
        return self.pixsizeEW

    # set pixsizeNS value
    def set_pixSizeNS(self, pixsizeNS):
        self.pixsizeNS = pixsizeNS
   
    # get pixsizeNS value
    def get_pixSizeNS(self):
        return self.pixsizeNS
        
    # set rows number
    def set_rows(self, rows):
        self.rows = rows
   
    # get rows number
    def get_rows(self):
        return self.rows

    # set cols number
    def set_cols(self, cols):
        self.cols = cols
   
    # get cols number
    def get_cols(self):
        return self.cols

    # set rotationA value
    def set_rotationA(self, rotationA):
        self.rotationA = rotationA
   
    # get rotationA value
    def get_rotationA(self):
        return self.rotationA
        
    # set rotationB value
    def set_rotationB(self, rotationB):
        self.rotationB = rotationB
   
    # get rotationB value
    def get_rotationB(self):
        return self.rotationB

    # check absence of axis rotations or pixel size differences
    def check_params(self, tolerance = 1e-06):
        
        # check if pixel size can be considered the same in the two axis directions
        if abs(abs(self.pixsizeEW) - abs(self.pixsizeNS))/abs(self.pixsizeNS) > tolerance :
            return False, 'Pixel sizes in x and y directions are different in raster'
            
        # check for the absence of axis rotations
        if abs(self.rotationA) > tolerance or abs(self.rotationB) > tolerance:
            return False, 'There should be no axis rotation in raster'
        
        return True, 'OK'

    # checks equivalence between the geographical parameters of two grids
    def geo_equiv(self, other, tolerance=1.0e-6): 

        if 2*(self.get_topLeftX() - other.get_topLeftX())/(self.get_topLeftX() + other.get_topLeftX()) > tolerance or \
            2*(self.get_topLeftY() - other.get_topLeftY())/(self.get_topLeftY() + other.get_topLeftY()) > tolerance or \
            2*(abs(self.get_pixSizeEW()) - abs(other.get_pixSizeEW()))/(abs(self.get_pixSizeEW()) + abs(other.get_pixSizeEW())) > tolerance or \
            2*(abs(self.get_pixSizeNS()) - abs(other.get_pixSizeNS()))/(abs(self.get_pixSizeNS()) + abs(other.get_pixSizeNS())) > tolerance or \
            self.get_rows() != other.get_rows() or self.get_cols() != other.get_cols() or \
            self.get_projection() != other.get_projection():
            return False
        else:
            return True


# exception for raster parameters
class Raster_Parameters_Errors(Exception):
    pass  


# exception for function input errors
class FunInp_Err(Exception):
    pass
    
    
# exception for output errors
class Output_Errors(Exception):
    pass
   
    
# Rectangular spatial domain class
class SpatialDomain(object):
    
    # class constructor
    def __init__(self, pt_init, pt_end): 
        self.pt_init = pt_init # lower-left corner of the domain, class: Point
        self.pt_end = pt_end # top-right corner of the domain, class: Point

    # get start point
    def get_start_point(self):
        return self.pt_init

    # get end point
    def get_end_point(self):
        return self.pt_end
    
    # get x range of spatial domain
    def get_xrange(self):
        return self.pt_end.x-self.pt_init.x

    # get y range of spatial domain
    def get_yrange(self):
        return self.pt_end.y-self.pt_init.y

    # get z range of spatial domain
    def get_zrange(self):
        return self.pt_end.z-self.pt_init.z

    # get horizontal area of spatial domain
    def get_horiz_area(self):
        return self.get_xrange()*self.get_yrange()


# check array input
def arr_check(curr_array):
    if len(curr_array.shape) != 3 and \
        curr_array.shape[2] != 2:
        return False
    return True


# Module principal class
class Grid(object):

    # alternative class constructor
    def __init__(self, grid_params = None, grid_data = None): 
                
        if grid_params is not None:
            # define lower-left corner as initial point
            pt_init = Point(grid_params.get_topLeftX(), grid_params.get_topLeftY() - abs(grid_params.get_pixSizeNS())*grid_params.get_rows())
            # define top-right corner as end point
            pt_end = Point(grid_params.get_topLeftX() + abs(grid_params.get_pixSizeEW())*grid_params.get_cols(), grid_params.get_topLeftY())
                    
            self.grid_domain = SpatialDomain(pt_init, pt_end)
        else:
            self.grid_domain = None
            
        self.grid_data = grid_data

    # set grid domain
    def set_grid_domain(self, pt_init, pt_end):
        self.grid_domain = SpatialDomain(pt_init, pt_end)

    # get grid domain
    def get_grid_domain(self):
        return self.grid_domain

    # set grid data
    def set_grid_data(self, data_array, level=None):
        try:
            if level is None:
                self.grid_data = data_array
            elif level == 0 or level == 1:
                self.grid_data[:,:, level] = data_array
            else:
                raise Raster_Parameters_Errors('Provided level for data array must be 0 or 1')
        except:
            raise Raster_Parameters_Errors('Unable to set grid data (function set_grid_data error)')

    # get grid data
    def get_grid_data(self, level=None):        
        try:
            if level is None:
                return self.grid_data
            elif level == 0 or level == 1:
                return self.grid_data[:, :, level]
            else:
                raise Raster_Parameters_Errors('Provided level for data array extraction must be 0 or 1')
        except:
            raise Raster_Parameters_Errors('Unable to extract grid data (function get_grid_data error)')

    # get row number of the grid domain
    def get_ylines_num(self):
        return np.shape(self.grid_data)[0]		

    # column number of the grid domain 		
    def get_xlines_num(self):
        return np.shape(self.grid_data)[1]		

    # returns the cell size of the gridded dataset in the x direction 
    def get_cellsize_x(self):
        return self.grid_domain.get_xrange()/float(self.get_xlines_num())

    # returns the cell size of the gridded dataset in the y direction 
    def get_cellsize_y(self):
        return self.grid_domain.get_yrange()/float(self.get_ylines_num())

    # returns the mean horizontal cell size 
    def get_cellsize_horiz_mean(self):
        return (self.get_cellsize_x()+self.get_cellsize_y())/2.0

    # converts from geographic to grid coordinates
    def geog2gridcoord(self, curr_Pt):
        currArrCoord_grid_i = (self.get_grid_domain().get_end_point().y - curr_Pt.y)/self.get_cellsize_y()
        currArrCoord_grid_j = (curr_Pt.x - self.get_grid_domain().get_start_point().x)/self.get_cellsize_x()
        return ArrCoord(currArrCoord_grid_i, currArrCoord_grid_j)

    # converts from grid to geographic coordinates
    def grid2geogcoord(self, currArrCoord):
        currPt_geogr_y = self.get_grid_domain().get_end_point().y - currArrCoord.i*self.get_cellsize_y()
        currPt_geogr_x = self.get_grid_domain().get_start_point().x + currArrCoord.j*self.get_cellsize_x()
        return Point(currPt_geogr_x, currPt_geogr_y)

    # calculates magnitude 
    def magnitude(self):

        if not arr_check(self.grid_data):
            raise FunInp_Err('input requires data array with three dimensions and 2-level third dimension')
            
        vx = self.grid_data[:, :, 0]
        vy = self.grid_data[:, :, 1]

        magnitude = np.sqrt(vx**2 + vy**2)

        magnitude_fld = Grid()
        magnitude_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        magnitude_fld.set_grid_data(magnitude)
        
        return magnitude_fld

    # calculates orientations 
    def orientations(self):        

        if not arr_check(self.grid_data):
            raise FunInp_Err('input requires data array with three dimensions and 2-level third dimension')
                     
        vx, vy = self.grid_data[:, :, 0], self.grid_data[:, :, 1]

        orientations = np.degrees(np.arctan2(vx, vy))
        orientations = np.where(orientations < 0.0, orientations + 360.0, orientations)
        
        orientations_fld = Grid()
        orientations_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        orientations_fld.set_grid_data(orientations)
        
        return orientations_fld

    # calculates divergence
    def divergence(self):

        if not arr_check(self.grid_data):
            raise FunInp_Err('input requires data array with three dimensions and 2-level third dimension')

        dvx_dx = np.gradient(self.grid_data[:, :, 0])[1]
        dvy_dy = -(np.gradient(self.grid_data[:, :, 1])[0])

        divergence_fld = Grid()
        divergence_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        divergence_fld.set_grid_data((dvx_dx + dvy_dy)/self.get_cellsize_horiz_mean())

        return divergence_fld

    # calculates curl module
    def curl_module(self):

        if not arr_check(self.grid_data):
            raise FunInp_Err('input requires data array with three dimensions and 2-level third dimension')
                     
        dvy_dx = np.gradient(self.grid_data[:, :, 1])[1]
        dvx_dy = -(np.gradient(self.grid_data[:, :, 0])[0])
        
        curl_fld = Grid()
        curl_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        curl_fld.set_grid_data((dvy_dx - dvx_dy)/self.get_cellsize_horiz_mean())
        
        return curl_fld

    # calculates magnitude gradient along x axis
    def grad_xaxis(self):	

        if not arr_check(self.grid_data):
            raise FunInp_Err('input requires data array with three dimensions and 2-level third dimension')
                      
        vx, vy = self.grid_data[:, :, 0], self.grid_data[:, :, 1]

        vect_magn = np.sqrt(vx**2 + vy**2)
        dm_dy, dm_dx = np.gradient(vect_magn)
        
        vect_xgrad_fld = Grid()
        vect_xgrad_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        vect_xgrad_fld.set_grid_data(dm_dx/self.get_cellsize_horiz_mean())
        
        return vect_xgrad_fld

    # calculates magnitude gradient along y axis
    def grad_yaxis(self):	

        if not arr_check(self.grid_data):
            raise FunInp_Err('input requires data array with three dimensions and 2-level third dimension')
                      
        vx, vy = self.grid_data[:, :, 0], self.grid_data[:, :, 1]

        vect_magn = np.sqrt(vx**2 + vy**2)
        dm_dy, dm_dx = np.gradient(vect_magn)
        dm_dy = - dm_dy
        
        vect_ygrad_fld = Grid()
        vect_ygrad_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        vect_ygrad_fld.set_grid_data(dm_dy/self.get_cellsize_horiz_mean())
        
        return vect_ygrad_fld

    # calculates gradient along flow lines
    def grad_flowlines(self):	

        if not arr_check(self.grid_data):
            raise FunInp_Err('input requires data array with three dimensions and 2-level third dimension')
                      
        vx, vy = self.grid_data[:, :, 0], self.grid_data[:, :, 1]
        
        dir_array = np.arctan2(vx, vy)
        
        vect_magn = np.sqrt(vx**2 + vy**2)
        dm_dy, dm_dx = np.gradient(vect_magn)
        dm_dy = - dm_dy
        
        velocity_gradient = dm_dx * np.sin(dir_array) + dm_dy * np.cos(dir_array)
        
        vect_magn_grad_fld = Grid()
        vect_magn_grad_fld.set_grid_domain(self.get_grid_domain().get_start_point(), self.get_grid_domain().get_end_point())
        vect_magn_grad_fld.set_grid_data(velocity_gradient/self.get_cellsize_horiz_mean())
        
        return vect_magn_grad_fld

    # returns the velocity components interpolated using the bilinear interpolation		
    def interpolate_level_bilinear(self, level, curr_Pt_gridcoord):
 
        try:
            v_00 = self.get_grid_data(level)[floor(curr_Pt_gridcoord.i - 0.5), floor(curr_Pt_gridcoord.j - 0.5)]
            v_10 = self.get_grid_data(level)[floor(curr_Pt_gridcoord.i - 0.5), ceil(curr_Pt_gridcoord.j - 0.5)]
            v_01 = self.get_grid_data(level)[ceil(curr_Pt_gridcoord.i - 0.5), floor(curr_Pt_gridcoord.j - 0.5)]    
            v_11 = self.get_grid_data(level)[ceil(curr_Pt_gridcoord.i - 0.5), ceil(curr_Pt_gridcoord.j - 0.5)]
        except:
            raise Raster_Parameters_Errors('Error in get_grid_data() function')
            
        delta_j_grid = curr_Pt_gridcoord.j - int(curr_Pt_gridcoord.j)
        
        if not (0 <= delta_j_grid < 1):
            raise Raster_Parameters_Errors('Delta j grid value must be between 0 and 1')
        
        if delta_j_grid >= 0.5:
            delta_j_grid = delta_j_grid - 0.5
        else:
            delta_j_grid = 0.5 + delta_j_grid             
        
        v_y0 = v_00 + (v_10-v_00)*delta_j_grid
        v_y1 = v_01 + (v_11-v_01)*delta_j_grid

        delta_i_grid = curr_Pt_gridcoord.i - int(curr_Pt_gridcoord.i)

        if not (0 <= delta_i_grid < 1):
            raise Raster_Parameters_Errors('Delta i grid value must be between 0 and 1')
        
        if delta_i_grid >= 0.5:
            delta_i_grid = delta_i_grid - 0.5
        else:
            delta_i_grid = 0.5 + delta_i_grid
            
        interp_v = v_y0 + (v_y1-v_y0)*delta_i_grid
                
        return interp_v  

    # check if a point pathline can be evaluated
    def include_point_location(self, currPoint): 
                    
        curr_Pt_gridcoord = self.geog2gridcoord(currPoint)
        
        if curr_Pt_gridcoord.i < 0.5 or curr_Pt_gridcoord.j < 0.5 or \
           curr_Pt_gridcoord.j > (self.get_xlines_num() - 0.5) or curr_Pt_gridcoord.i > (self.get_ylines_num() - 0.5):
            return False
 
        if len(self.grid_data.shape) != 3 and \
           self.grid_data.shape[2] != 2:
            raise FunInp_Err('input requires data array with three dimensions and 2-level third dimension')
                    
        for n in range(2):
            if np.isnan(self.get_grid_data(n)[floor(curr_Pt_gridcoord.i - 0.5), floor(curr_Pt_gridcoord.j - 0.5)]) or \
               np.isnan(self.get_grid_data(n)[ceil(curr_Pt_gridcoord.i - 0.5), floor(curr_Pt_gridcoord.j - 0.5)]) or \
               np.isnan(self.get_grid_data(n)[floor(curr_Pt_gridcoord.i - 0.5), ceil(curr_Pt_gridcoord.j - 0.5)]) or \
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
        
        K2_Pt = Point(curr_Pt.x + (0.25)*delta_time*K1_vx, curr_Pt.y + (0.25)*delta_time*K1_vy) 
        K2_vx, K2_vy = self.point_velocity(K2_Pt)
        if K2_vx is None or K2_vy is None:
            return None, None

        K3_Pt = Point(curr_Pt.x + (3.0/32.0)*delta_time*K1_vx + (9.0/32.0)*delta_time*K2_vx, \
                      curr_Pt.y + (3.0/32.0)*delta_time*K1_vy + (9.0/32.0)*delta_time*K2_vy) 
        K3_vx, K3_vy = self.point_velocity(K3_Pt)
        if K3_vx is None or K3_vy is None:
            return None, None

        K4_Pt = Point(curr_Pt.x + (1932.0/2197.0)*delta_time*K1_vx - (7200.0/2197.0)*delta_time*K2_vx + (7296.0/2197.0)*delta_time*K3_vx, \
                      curr_Pt.y + (1932.0/2197.0)*delta_time*K1_vy - (7200.0/2197.0)*delta_time*K2_vy + (7296.0/2197.0)*delta_time*K3_vy) 
        K4_vx, K4_vy = self.point_velocity(K4_Pt)
        if K4_vx is None or K4_vy is None:
            return None, None
        
        K5_Pt = Point(curr_Pt.x + (439.0/216.0)*delta_time*K1_vx - (8.0)*delta_time*K2_vx + (3680.0/513.0)*delta_time*K3_vx - (845.0/4104.0)*delta_time*K4_vx, \
                      curr_Pt.y + (439.0/216.0)*delta_time*K1_vy - (8.0)*delta_time*K2_vy + (3680.0/513.0)*delta_time*K3_vy - (845.0/4104.0)*delta_time*K4_vy) 
        K5_vx, K5_vy = self.point_velocity(K5_Pt)
        if K5_vx is None or K5_vy is None:
            return None, None
            
        K6_Pt = Point(curr_Pt.x - (8.0/27.0)*delta_time*K1_vx + (2.0)*delta_time*K2_vx - (3544.0/2565.0)*delta_time*K3_vx + (1859.0/4104.0)*delta_time*K4_vx - (11.0/40.0)*delta_time*K5_vx, \
                      curr_Pt.y - (8.0/27.0)*delta_time*K1_vy + (2.0)*delta_time*K2_vy - (3544.0/2565.0)*delta_time*K3_vy + (1859.0/4104.0)*delta_time*K4_vy - (11.0/40.0)*delta_time*K5_vy) 
        K6_vx, K6_vy = self.point_velocity(K6_Pt)
        if K6_vx is None or K6_vy is None:
            return None, None            
         
        rkf_4o_x = curr_Pt.x + delta_time*((25.0/216.0)*K1_vx + (1408.0/2565.0)*K3_vx + (2197.0/4104.0)*K4_vx - (1.0/5.0)*K5_vx)
        rkf_4o_y = curr_Pt.y + delta_time*((25.0/216.0)*K1_vy + (1408.0/2565.0)*K3_vy + (2197.0/4104.0)*K4_vy - (1.0/5.0)*K5_vy)
        temp_Pt = Point(rkf_4o_x, rkf_4o_y)

        interp_Pt_x = curr_Pt.x + delta_time*((16.0/135.0)*K1_vx + (6656.0/12825.0)*K3_vx + (28561.0/56430.0)*K4_vx - (9.0/50.0)*K5_vx + (2.0/55.0)*K6_vx)
        interp_Pt_y = curr_Pt.y + delta_time*((16.0/135.0)*K1_vy + (6656.0/12825.0)*K3_vy + (28561.0/56430.0)*K4_vy - (9.0/50.0)*K5_vy + (2.0/55.0)*K6_vy)
        interp_Pt = Point(interp_Pt_x, interp_Pt_y)
        
        interp_PT_error_estimate = interp_Pt.distance(temp_Pt)
        
        return interp_Pt, interp_PT_error_estimate

    # writes ESRI ascii grid
    def write_esrigrid(self, outgrid_fn, esri_nullvalue=-99999, level=0):

        outgrid_fn = str(outgrid_fn)
        
        # checking existence of output slope grid
        if os.path.exists(outgrid_fn):
            raise Output_Errors("Output grid '%s' already exists" % outgrid_fn)

        try:
            outputgrid = open(outgrid_fn, 'w')  # create the output ascii file
        except:
            raise Output_Errors("Unable to create output grid '%s'" % outgrid_fn)
       
        if outputgrid is None:
            raise Output_Errors("Unable to create output grid '%s'" % outgrid_fn)

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
                            outputgrid.write('%.8f ' % (esrigrid_outvalues[i,j,level]))
                    outputgrid.write('\n')
        elif len(self.grid_data.shape) == 2:
            for i in range(0, self.get_ylines_num()):
                    for j in range(0, self.get_xlines_num()):
                            outputgrid.write('%.8f ' % (esrigrid_outvalues[i,j]))
                    outputgrid.write('\n')
                    
        outputgrid.close()

        return True        


