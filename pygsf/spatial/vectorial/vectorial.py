# -*- coding: utf-8 -*-


from ...defaults.orientations import *
from ...mathematics.vectors import *


class Point(object):
    """
    Cartesian point.
    Dimensions: 3D (space)
    """

    def __init__(self, x: [int, float], y: [int, float], z: [int, float]):
        """
        Construct a Point instance.
        """

        vals = [x, y, z]
        if any(map(lambda val: not isinstance(val, (int, float)), vals)):
            raise VectorInputException("Input values must be integer of float")
        elif not all(map(isfinite, vals)):
            raise VectorInputException("Input values must be finite")
        else:
            self._a = array(vals, dtype=np.float64)

    def __repr__(self) -> str:

        return "Point({:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z)

    def __eq__(self, another: 'Point') -> Optional[bool]:
        """
        Return True if objects are equal.

        Example:
          >>> Point(1., 1., 1.) == Point(1, 1, 1)
          True
          >>> Point(1., 1., 1.) == Point(1, 1, -1)
          False
        """

        if not isinstance(another, Point):
            raise VectorInputException("Variables must be of the same type")
        else:
            return all([
                self.x == another.x,
                self.y == another.y,
                self.z == another.z])

    def __ne__(self, another: 'Point') -> Optional[bool]:
        """
        Return False if objects are equal.

        Example:
          >>> Point(1., 1., 1.) != Point(0., 0., 0.)
          True
        """

        if not isinstance(another, Point):
            return None
        else:
            return not (self == another)

    @property
    def a(self) -> 'np.array':
        """
        Return a copy of the object inner array.

        :return: double array of x, y, z values

        Examples:
          >>> Point(4, 3, 7).a
          array([4., 3., 7.])
        """

        return np.copy(self._a)

    @property
    def x(self) -> float:
        """
        Return x value

        Example:
          >>> Point(1.5, 1, 1).x
          1.5
        """

        return self.a[0]

    @property
    def y(self) -> float:
        """
        Return y value

        Example:
          >>> Point(1.5, 3.0, 1).y
          3.0
        """
        return self.a[1]

    @property
    def z(self) -> float:
        """
        Return z value

        Example:
          >>> Point(1.5, 3.2, 41.).z
          41.0
        """
        return self.a[2]

    def toXYZ(self) -> Tuple[float, float, float]:
        """
        Returns the spatial components as a tuple of three values.

        :return: the spatial components (x, y, z).
        :rtype: a tuple of three floats.

        Examples:
          >>> Point(1, 0, 3).toXYZ()
          (1.0, 0.0, 3.0)
        """

        return self.x, self.y, self.z

    def toArray(self) -> 'array':
        """
        Return a double Numpy array representing the point values.

        :return: Numpy array

        Examples:
          >>> Point(1, 2, 3).toArray()
          array([1., 2., 3.])
        """

        return self.a

    def pXY(self) -> 'Point':
        """
        Projection on the x-y plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pXY()
          Point(2.0000, 3.0000, 0.0000)
        """

        return self.__class__(self.x, self.y, 0.0)

    def pXZ(self) -> 'Point':
        """
        Projection on the x-z plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pXZ()
          Point(2.0000, 0.0000, 4.0000)
        """

        return self.__class__(self.x, 0.0, self.z)

    def pYZ(self) -> 'Point':
        """
        Projection on the y-z plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pYZ()
          Point(0.0000, 3.0000, 4.0000)
        """

        return self.__class__(0.0, self.y, self.z)

    @property
    def len3D(self) -> float:
        """
        Spatial distance of the point from the axis origin.

        :return: distance
        :rtype: float

        Examples:
          >>> Point(4.0, 3.0, 0.0).len3D
          5.0
        """

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    @property
    def len2D(self) -> float:
        """
        2D distance of the point from the axis origin.

        Example:
          >>> Point(3, 4, 0).len2D
          5.0
          >>> Point(12, 5, 3).len2D
          13.0
        """

        return sqrt(self.x * self.x + self.y * self.y)

    def deltaX(self, another: 'Point') -> Optional[float]:
        """
        Delta between x components of two Point Instances.

        :return: difference value
        :rtype: float

        Examples:
          >>> Point(1, 2, 3).deltaX(Point(4, 7, 1))
          3.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return another.x - self.x

    def deltaY(self, another: 'Point') -> Optional[float]:
        """
        Delta between y components of two Point Instances.

        :return: difference value
        :rtype: float

        Examples:
          >>> Point(1, 2, 3).deltaY(Point(4, 7, 1))
          5.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return another.y - self.y

    def deltaZ(self, another: 'Point') -> Optional[float]:
        """
        Delta between x components of two Point Instances.

        :return: difference value
        :rtype: float

        Examples:
          >>> Point(1, 2, 3).deltaZ(Point(4, 7, 1))
          -2.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return another.z - self.z

    def dist3DWith(self, another: 'Point') -> Optional[float]:
        """
        Calculate Euclidean spatial distance between two points.

        Examples:
          >>> Point(1., 1., 1.).dist3DWith(Point(4., 5., 1,))
          5.0
          >>> Point(1, 1, 1).dist3DWith(Point(4, 5, 1))
          5.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2 + (self.z - another.z) ** 2)

    def dist2DWith(self, another: 'Point') -> Optional[float]:
        """
        Calculate horizontal (2D) distance between two points.

        Examples:
          >>> Point(1., 1., 1.).dist2DWith(Point(4., 5., 7.))
          5.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def scale(self, scale_factor: [int, float]) -> Optional['Point']:
        """
        Create a scaled object.

        Example;
          >>> Point(1, 0, 1).scale(2.5)
          Point(2.5000, 0.0000, 2.5000)
          >>> Point(1, 0, 1).scale(np.nan) is None
          True
          >>> Point(1, 0, 1).scale(np.inf) is None
          True
        """

        if not isinstance(scale_factor, (int, float)):
            return None
        elif not isfinite(scale_factor):
            return None
        else:
            x, y, z = arr2tuple(self.a * scale_factor)
            return self.__class__(x, y, z)

    def invert(self) -> 'Point':
        """
        Create a new object with inverted direction.

        Examples:
          >>> Point(1, 1, 1).invert()
          Point(-1.0000, -1.0000, -1.0000)
          >>> Point(2, -1, 4).invert()
          Point(-2.0000, 1.0000, -4.0000)
        """

        return self.scale(-1)

    def isCoinc(self, another: 'Point', tolerance: [int, float] = MIN_SEPARATION_THRESHOLD) -> Optional[bool]:
        """
        Check spatial coincidence of two points

        Example:
          >>> Point(1., 0., -1.).isCoinc(Point(1., 1.5, -1.))
          False
          >>> Point(1., 0., 0.).isCoinc(Point(1., 0., 0.))
          True
          >>> Point(1.2, 7.4, 1.4).isCoinc(Point(1.2, 7.4, 1.4))
          True
          >>> Point(1.2, 7.4, 1.4).isCoinc(Point(1.2, 7.4, 1.4), tolerance=np.nan) is None
          True
        """

        if not isinstance(another, Point):
            return None
        elif not isinstance(tolerance, (int, float)):
            return None
        elif not isfinite(tolerance):
            return None
        else:
            distance_2d = self.dist2DWith(another)
            if np.isnan(distance_2d) or distance_2d > tolerance:
                return False
            else:
                distance_3d = self.dist3DWith(another)
                if np.isnan(distance_3d) or distance_3d > tolerance:
                    return False
                else:
                    return True

    def shift(self, sx: float, sy: float, sz: float) -> Optional['Point']:
        """
        Create a new object shifted by given amount from the self instance.

        Example:
          >>> Point(1, 1, 1).shift(0.5, 1., 1.5)
          Point(1.5000, 2.0000, 2.5000)
          >>> Point(1, 2, -1).shift(0.5, 1., 1.5)
          Point(1.5000, 3.0000, 0.5000)
          >>> Point(1, 2, -1).shift(0.5, np.nan, 1.5) is None
          True
       """

        vals = [sx, sy, sz]
        if not all(map(lambda val: isinstance(val, (int, float)), vals)):
            return None
        elif not all(map(isfinite, vals)):
            return None
        else:
            return self.__class__(self.x + sx, self.y + sy, self.z + sz)

    def asVect(self) -> 'Vect':
        """
        Create a vector based on the point coordinates

        Example:
          >>> Point(1, 1, 0).asVect()
          Vect(1.0000, 1.0000, 0.0000)
          >>> Point(0.2, 1, 6).asVect()
          Vect(0.2000, 1.0000, 6.0000)
        """

        return Vect(self.x, self.y, self.z)


class CPlane(object):
    """
    Cartesian plane.
    Expressed by equation:
    ax + by + cz + d = 0

    Note: Plane is locational - its position in space is defined.
    This contrast with PPlane, defined just by its attitude, but with undefined position

    """

    def __init__(self, a, b, c, d):

        self._a = float(a)
        self._b = float(b)
        self._c = float(c)
        self._d = float(d)

    @property
    def a(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> CPlane(1, 0, 0, 2).a
          1.0
        """

        return self._a

    @property
    def b(self):
        """
        Return b coefficient of a Plane instance.

        Example:
          >>> CPlane(1, 4, 0, 2).b
          4.0
        """

        return self._b

    @property
    def c(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> CPlane(1, 0, 5.4, 2).c
          5.4
        """

        return self._c

    @property
    def d(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> CPlane(1, 0, 0, 2).d
          2.0
        """

        return self._d

    @property
    def v(self):
        """
        Return coefficients of a Plane instance.

        Example:
          >>> CPlane(1, 1, 7, -4).v
          (1.0, 1.0, 7.0, -4.0)
        """
        return self.a, self.b, self.c, self.d

    @classmethod
    def fromPoints(cls, pt1, pt2, pt3):
        """
        Create a Plane from three given Point instances.

        Example:
          >>> CPlane.fromPoints(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0))
          CPlane(0.0000, 0.0000, 1.0000, 0.0000)
          >>> CPlane.fromPoints(Point(0, 0, 0), Point(0, 1, 0), Point(0, 0, 1))
          CPlane(1.0000, 0.0000, 0.0000, 0.0000)
        """

        matr_a = array(
            [[pt1.y, pt1.z, 1],
            [pt2.y, pt2.z, 1],
            [pt3.y, pt3.z, 1]])

        matr_b = - array(
            [[pt1.x, pt1.z, 1],
            [pt2.x, pt2.z, 1],
            [pt3.x, pt3.z, 1]])

        matr_c = array(
            [[pt1.x, pt1.y, 1],
            [pt2.x, pt2.y, 1],
            [pt3.x, pt3.y, 1]])

        matr_d = - array(
            [[pt1.x, pt1.y, pt1.z],
            [pt2.x, pt2.y, pt2.z],
            [pt3.x, pt3.y, pt3.z]])

        return cls(
            np.linalg.det(matr_a),
            np.linalg.det(matr_b),
            np.linalg.det(matr_c),
            np.linalg.det(matr_d))

    def __repr__(self):

        return "CPlane({:.4f}, {:.4f}, {:.4f}, {:.4f})".format(*self.v)

    def normVersor(self):
        """
        Return the versor normal to the cartesian plane.

        Examples:
          >>> CPlane(0, 0, 5, -2).normVersor()
          Vect(0.0000, 0.0000, 1.0000)
          >>> CPlane(0, 7, 0, 5).normVersor()
          Vect(0.0000, 1.0000, 0.0000)
        """

        return Vect(self.a, self.b, self.c).versor()

    def toPoint(self):
        """
        Returns a point lying in the plane (non-unique solution).

        Examples:
          >>> CPlane(0, 0, 1, -1).toPoint()
          Point(0.0000, 0.0000, 1.0000)
        """

        point = Point(*point_solution(array([[self.a, self.b, self.c]]),
                                      array([-self.d])))
        return point

    def intersVersor(self, another):
        """
        Return intersection versor for two intersecting planes.

        >>> a = CPlane(1, 0, 0, 0)
        >>> b = CPlane(0, 0, 1, 0)
        >>> a.intersVersor(b)
        Vect(0.0000, -1.0000, 0.0000)
        """

        return self.normVersor().vCross(another.normVersor()).versor()

    def intersPoint(self, another):
        """
        Return point on intersection line (obviously non-unique solution)
        for two planes.

        >>> a = CPlane(1, 0, 0, 0)
        >>> b = CPlane(0, 0, 1, 0)
        >>> a.intersPoint(b)
        Point(0.0000, 0.0000, 0.0000)
        """

        # find a point lying on the intersection line (this is a non-unique solution)
        a = array([[self.a, self.b, self.c], [another.a, another.b, another.c]])
        b = array([-self.d, -another.d])
        x, y, z = point_solution(a, b)

        return Point(x, y, z)

    def isPointInPlane(self, pt):
        """
          Check whether a point lie in a plane.

          >>> pl = CPlane(0, 0, 1, 0)
          >>> pt = Point(0, 1, 0)
          >>> pl.isPointInPlane(pt)
          True
        """

        if abs(self.a * pt.x + self.b * pt.y + self.c * pt.z + self.d) < MIN_SCALAR_VALUE:
            return True
        else:
            return False

    def angle(self, another):
        """
        Calculate angle (in degrees) between two planes.

        Examples:
          >>> CPlane(1,0,0,0).angle(CPlane(0,1,0,0))
          90.0
          >>> CPlane(1,0,0,0).angle(CPlane(1,0,1,0))
          45.0
          >>> CPlane(1,0,0,0).angle(CPlane(1,0,0,0))
          0.0
        """

        angle_degr = self.normVersor().angle(another.normVersor())
        if abs(angle_degr) < MIN_ANGLE_DEGR_VALUE:
            angle_degr = 0.0
        elif angle_degr > 90.0:
            angle_degr = 180.0 - angle_degr

        return angle_degr

    def isSubParallel(self, another, angle_tolerance=PLANE_ANGLE_THRESHOLD):
        """
        Check that two Plane are sub-parallel

        :param another: a Plane instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

         Examples:
          >>> CPlane(1,0,0,0).isSubParallel(CPlane(1,0,0,0))
          True
          >>> CPlane(1,0,0,0).isSubParallel(CPlane(1,0,1,0))
          False
        """

        return self.angle(another) < angle_tolerance


# TODO class Segment
"""
Defined by a couple of points
"""

# TODO class Profile
"""
A 1D position values plus M(s) values
"""

# TODO class Path
"""
The trajectory of a Point with time
"""


if __name__ == "__main__":

    import doctest
    doctest.testmod()
