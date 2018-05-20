# -*- coding: utf-8 -*-


from .scalars import *


def arrToTuple(arr1D: 'array[Number]') -> Tuple[float, ...]:
    """
    Modified from: https://stackoverflow.com/questions/10016352/convert-numpy-array-to-tuple
    Works just for 1D arrays

    :param arr1D: the 1D-arrays whose components have to be extracted
    :type arr1D: numpy array
    :return: a tuple derived from the array values extraction
    :rtype: tuple of float

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
    :rtype: list of Floats

    Examples:
      >>> toFloats([1, 2, 3])
      [1.0, 2.0, 3.0]
    """

    return [float(item) for item in iterable_obj]


def arraysAreClose(a_array: 'array[Number]', b_array: 'array[Number]',
                   rtol: float=1e-012, atol: float=1e-12, equal_nan: bool=False, equal_inf: bool=False) -> bool:
    """
    Check for equivalence between two numpy arrays.

    :param a_array: first array to be compared
    :type a_array: numpy array
    :param b_array: second array to be compared with the first one
    :type b_array: numpy array
    :param rtol: relative tolerance
    :type rtol:
    :param atol: absolute tolerance
    :type atol:
    :param equal_nan: consider nan values equivalent or not
    :type equal_nan:
    :param equal_inf: consider inf values equivalent or not
    :type equal_inf:
    :return: whether the two arrays are close as component values
    :rtype: bool

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


def pointSolution(a_array: 'array[Number]', b_array: 'array[Number]'):
    """
    Finds a non-unique solution for a set of linear equations.

    :param a_array:
    :type a_array: numpy array
    :param b_array:
    :type b_array: numpy array
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
    :type xyz_array: numpy array
    :return:
    :rtype:

    Examples:
    """

    try:
        result = np.linalg.svd(xyz_array)
    except:
        result = None

    return dict(result=result)


def divergence_2D(fld_x: 'array', fld_y: 'array', cell_size_x: Number, cell_size_y: Number) -> 'array':
    """
    Calculates the divergence from two 2D arrays:
    the first represents the vector field x component, the second the vector field y component.

    :param fld_x: vector field x component
    :type fld_x: np.array
    :param fld_y: vector field y component
    :type fld_y: np.array
    :return: divergence field
    :rtype: np.array

    Examples:
    """

    grad_x = np.gradient(fld_x, edge_order=2, axis=1)
    grad_y = np.gradient(fld_y, edge_order=2, axis=0)

    dfx_dx = grad_x / cell_size_x
    dfy_dy = grad_y / cell_size_y

    return dfx_dx - dfy_dy


if __name__ == "__main__":

    import doctest
    doctest.testmod()
