# -*- coding: utf-8 -*-


from ..defaults.mathematics import *
from ..defaults.typing import *

from .scalars import *


def arrToTuple(arr1D: 'array[Numbers]') -> Tuple[float, ...]:
    """
    Modified from: https://stackoverflow.com/questions/10016352/convert-numpy-array-to-tuple
    Works just for 1D arrays

    :param arr1D: the 1D-arrays whose components have to be extracted
    :type arr1D: numpy array
    :return: a tuple derived from the array values extraction
    :rtype: tuple of float

    Examples:
      >>> arr = array([1,2,3,4,5])
      >>> arrToTuple(arr)
      (1.0, 2.0, 3.0, 4.0, 5.0)
    """

    return tuple(map(float, arr1D))


def toFloats(iterable_obj: Sequence[Numbers]) -> List[float]:
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


def arraysAreClose(a_array: 'array[Numbers]', b_array: 'array[Numbers]',
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

    are_equal = []
    for a, b in np.nditer([a_array, b_array]):
        are_equal.append(areClose(a.item(0), b.item(0), rtol=rtol, atol=atol, equal_nan=equal_nan, equal_inf=equal_inf))

    return all(are_equal)


def pointSolution(a_array: 'array[Numbers]', b_array: 'array[Numbers]'):
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


if __name__ == "__main__":

    import doctest
    doctest.testmod()
