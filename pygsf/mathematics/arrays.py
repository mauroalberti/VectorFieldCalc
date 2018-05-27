# -*- coding: utf-8 -*-


from .scalars import *
from .exceptions import InputValuesException


def arrToTuple(arr1D: 'array[Number]') -> Tuple[float, ...]:
    """
    Modified from: https://stackoverflow.com/questions/10016352/convert-numpy-array-to-tuple
    Works just for 1D arrays

    :param arr1D: the 1D-arrays whose components have to be extracted.
    :type arr1D: numpy array.
    :return: a tuple derived from the array values extraction.
    :rtype: tuple of float.

    Examples:
      >>> levels = array([1,2,3,4,5])
      >>> arrToTuple(levels)
      (1.0, 2.0, 3.0, 4.0, 5.0)
    """

    return tuple(map(float, arr1D))


def toFloats(iterable_obj: Sequence[Number]) -> List[float]:
    """
    Converts an iterable object storing float-compatible values to a list of floats.

    :param iterable_obj:
    :type iterable_obj:
    :return:
    :rtype: list of Floats.

    Examples:
      >>> toFloats([1, 2, 3])
      [1.0, 2.0, 3.0]
    """

    return [float(item) for item in iterable_obj]


def arraysAreClose(a_array: 'array[Number]', b_array: 'array[Number]',
                   rtol: float=1e-012, atol: float=1e-12, equal_nan: bool=False, equal_inf: bool=False) -> bool:
    """
    Check for equivalence between two numpy arrays.

    :param a_array: first array to be compared.
    :type a_array: numpy array.
    :param b_array: second array to be compared with the first one.
    :type b_array: numpy array.
    :param rtol: relative tolerance.
    :type rtol:
    :param atol: absolute tolerance.
    :type atol:
    :param equal_nan: consider nan values equivalent or not.
    :type equal_nan:
    :param equal_inf: consider inf values equivalent or not.
    :type equal_inf:
    :return: whether the two arrays are close as component values.
    :rtype: bool.

    Examples:
      >>> arraysAreClose(array([1,2,3]), array([1,2,3]))
      True
      >>> arraysAreClose(array([[1,2,3], [4, 5, 6]]), array([1,2,3]))
      False
      >>> arraysAreClose(array([[1,2,3], [4,5,6]]), array([[1,2,3], [4,5,6]]))
      True
      >>> arraysAreClose(array([[1,2,np.nan], [4,5,6]]), array([[1,2,np.nan], [4,5,6]]))
      False
      >>> arraysAreClose(array([[1,2,np.nan], [4,5,6]]), array([[1,2,np.nan], [4,5,6]]), equal_nan=True)
      True
    """
    if a_array.shape != b_array.shape:
        return False

    are_close = []
    for a, b in np.nditer([a_array, b_array]):
        are_close.append(areClose(a.item(0), b.item(0), rtol=rtol, atol=atol, equal_nan=equal_nan, equal_inf=equal_inf))

    return all(are_close)


def interp_bilinear(arr: 'array', i: Number, j: Number) -> Optional[float]:
    """
    Interpolate the z value at a point, given its array coordinates.
    Interpolation method: bilinear.

    0,0   0,1

    1,0,  1,1

    :param arr: array with values for which the interpolation will be made.
    :type arr: Numpy array.
    :param i: i array index of the point (may be fractional).
    :type i: Number.
    :param j: j array index of the point (may be fractional).
    :type j: Number.
    :return: interpolated z value.
    :rtype: optional float.
    """

    if i < 0.5 or j < 0.5:
        return None

    i_max, j_max = arr.shape

    if i > i_max - 0.5 or j > j_max - 0.5:
        return None

    loc_cellcent_i = i - 0.5
    loc_cellcent_j = j - 0.5

    grid_val_00 = arr[int(floor(loc_cellcent_i)), int(floor(loc_cellcent_j))]
    grid_val_01 = arr[int(floor(loc_cellcent_i)), int(ceil(loc_cellcent_j))]
    grid_val_10 = arr[int(ceil(loc_cellcent_i)), int(floor(loc_cellcent_j))]
    grid_val_11 = arr[int(ceil(loc_cellcent_i)), int(ceil(loc_cellcent_j))]

    delta_i = loc_cellcent_i - floor(loc_cellcent_i)
    delta_j = loc_cellcent_j - floor(loc_cellcent_j)

    grid_val_y0 = grid_val_00 + (grid_val_10 - grid_val_00) * delta_i
    grid_val_y1 = grid_val_01 + (grid_val_11 - grid_val_01) * delta_i

    grid_val_interp = grid_val_y0 + (grid_val_y1 - grid_val_y0) * delta_j

    return grid_val_interp


def gradient_x(fld: 'array', cell_size_x: Number) -> 'array':
    """
    Calculates the array gradient along the x axis.

    :param fld: array.
    :type fld: np.array.
    :param cell_size_x: the cell spacing in the x direction.
    :type cell_size_x: Number.
    :return: gradient field.
    :rtype: np.array.

    Examples:
    """

    return np.gradient(fld, edge_order=2, axis=1) / cell_size_x


def gradient_y(fld: 'array', cell_size_y: Number) -> 'array':
    """
    Calculates the array gradient along the y axis.

    :param fld: array.
    :type fld: np.array.
    :param cell_size_y: the cell spacing in the y direction.
    :type cell_size_y: Number.
    :return: gradient field.
    :rtype: np.array.

    Examples:
    """

    return np.gradient(fld, edge_order=2, axis=0) / cell_size_y


def magnitude_2D(fld_x: 'array', fld_y: 'array') -> 'array':
    """
    Calculates the magnitude given two 2D arrays:
    the first represents the vector field x component, the second the vector field y component.

    :param fld_x: vector field x component.
    :type fld_x: np.array.
    :param fld_y: vector field y component.
    :type fld_y: np.array.
    :return: magnitude field.
    :rtype: np.array.

    Examples:
    """

    return np.sqrt(fld_x ** 2 + fld_y ** 2)


def magnitude_gradient(fld_x: 'array', fld_y: 'array', dir_cell_sizes: List[Number], axis: str='') -> List['array']:
    """
    Calculates the magnitude gradient along the x-direction, given two 2D arrays:
    the first represents the vector field x component, the second the vector field y component.

    :param fld_x: vector field x component.
    :type fld_x: np.array.
    :param fld_y: vector field y component.
    :type fld_y: np.array.
    :param dir_cell_sizes: list of cell spacing(s) in the considered direction(s).
    :type dir_cell_sizes: list of Number(s).
    :param axis: declares the axis ('x' or 'y') or the axes('', i.e., empty string) for both x and y directions.
    :type axis: str.
    :return: magnitude gradient field(s) along the considered direction.
    :rtype: list of np.array.
    :raises: InputValuesException.

    Examples:
    """

    magn = magnitude_2D(fld_x, fld_y)
    if axis == 'x':
        return [gradient_x(magn, dir_cell_sizes[0])]
    elif axis == 'y':
        return [gradient_y(magn, dir_cell_sizes[0])]
    elif axis == '':
        return [gradient_x(magn, dir_cell_sizes[0]), gradient_y(magn, dir_cell_sizes[1])]
    else:
        raise InputValuesException("Axis must be 'x' or 'y' or '' (for both x and y). '{}' given".format(axis))


def orientations_rad(fld_x: 'array', fld_y: 'array') -> 'array':
    """
    Calculates the orientations (as radians) given two 2D arrays:
    the first represents the vector field x component, the second the vector field y component.

    :param fld_x: vector field x component.
    :type fld_x: np.array.
    :param fld_y: vector field y component.
    :type fld_y: np.array.
    :return: orientation field, in radians.
    :rtype: np.array.

    Examples:
    """

    azimuth_rad = np.arctan2(fld_x, fld_y)
    azimuth_rad = np.where(azimuth_rad < 0.0, azimuth_rad + np.pi, azimuth_rad)

    return azimuth_rad


def orientations_degr(fld_x: 'array', fld_y: 'array') -> 'array':
    """
    Calculates the orientations (as decimal degrees) given two 2D arrays:
    the first represents the vector field x component, the second the vector field y component.

    :param fld_x: vector field x component.
    :type fld_x: np.array.
    :param fld_y: vector field y component.
    :type fld_y: np.array.
    :return: orientation field, in decimal degrees.
    :rtype: np.array.

    Examples:
    """

    return np.degrees(orientations_rad(fld_x, fld_y))


def divergence_2D(fld_x: 'array', fld_y: 'array', cell_size_x: Number, cell_size_y: Number) -> 'array':
    """
    Calculates the divergence from two 2D arrays:
    the first represents the vector field x component, the second the vector field y component.

    :param fld_x: vector field x component.
    :type fld_x: np.array.
    :param fld_y: vector field y component.
    :type fld_y: np.array.
    :param cell_size_x: the cell spacing in the x direction.
    :type cell_size_x: Number.
    :param cell_size_y: the cell spacing in the y direction.
    :type cell_size_y: Number.
    :return: divergence field.
    :rtype: np.array.

    Examples:
    """

    dfx_dx = gradient_x(fld_x, cell_size_x)
    dfy_dy = gradient_y(fld_y, cell_size_y)

    return dfx_dx - dfy_dy


def curl_mod(fld_x: 'array', fld_y: 'array', cell_size_x: Number, cell_size_y: Number) -> 'array':
    """
    Calculates the curl from two 2D arrays:
    the first represents the vector field x component, the second the vector field y component.

    :param fld_x: vector field x component.
    :type fld_x: np.array.
    :param fld_y: vector field y component.
    :type fld_y: np.array.
    :param cell_size_x: the cell spacing in the x direction.
    :type cell_size_x: Number.
    :param cell_size_y: the cell spacing in the y direction.
    :type cell_size_y: Number.
    :return: curl field.
    :rtype: np.array.

    Examples:
    """

    dfy_dx = gradient_x(fld_y, cell_size_x)
    dfx_dy = gradient_y(fld_x, cell_size_y)

    return dfy_dx - dfx_dy


def gradient_flowlines(fld_x: 'array', fld_y: 'array', cell_size_x: Number, cell_size_y: Number) -> 'array':
    """
    Calculates gradient along flow lines.

    :param fld_x: vector field x component.
    :type fld_x: np.array.
    :param fld_y: vector field y component.
    :type fld_y: np.array.
    :param cell_size_x: the cell spacing in the x direction.
    :type cell_size_x: Number.
    :param cell_size_y: the cell spacing in the y direction.
    :type cell_size_y: Number.
    :return: the flowline gradient field
    :rtype: np.array.
    """

    orien_rad = orientations_rad(fld_x, fld_y)

    dm_dx, dm_dy = magnitude_gradient(
        fld_x=fld_x,
        fld_y=fld_y,
        dir_cell_sizes=[cell_size_x, cell_size_y])
    dm_dy = - dm_dy

    velocity_gradient = dm_dx * np.sin(orien_rad) + dm_dy * np.cos(orien_rad)

    return velocity_gradient


def pointSolution(a_array: 'array[Number]', b_array: 'array[Number]'):
    """
    Finds a non-unique solution for a set of linear equations.

    :param a_array:
    :type a_array: numpy array.
    :param b_array:
    :type b_array: numpy array.
    :return:
    :rtype:

    Examples:
    """

    try:
        return np.linalg.lstsq(a_array, b_array, rcond=None)[0]
    except:
        return None, None, None


def xyzSvd(xyz_array) -> dict:
    """
    Calculates the SVD solution given a Numpy array.

    # modified after:
    # http://stackoverflow.com/questions/15959411/best-fit-plane-algorithms-why-different-results-solved

    :param xyz_array:
    :type xyz_array: numpy array.
    :return:
    :rtype:

    Examples:
    """

    try:
        result = np.linalg.svd(xyz_array)
    except:
        result = None

    return dict(result=result)


if __name__ == "__main__":

    import doctest
    doctest.testmod()
