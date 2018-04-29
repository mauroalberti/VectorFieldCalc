# -*- coding: utf-8 -*-


import numpy as np

from .scalars import *


array = np.array


def arr2tuple(arr1D) -> tuple:
    """
    Modified from: https://stackoverflow.com/questions/10016352/convert-numpy-array-to-tuple
    Works just for 1D arrays

    :param arr1D:
    :return: tuple of float values

    Examples:
      >>> arr = array([1,2,3,4,5])
      >>> arr2tuple(arr)
      (1.0, 2.0, 3.0, 4.0, 5.0)
    """

    return tuple(map(float, arr1D))


def to_floats(iterable_obj):
    """
    Converts an iterable object storing float-compatible values to a list of floats.

    :param iterable_obj:
    :return: list of Floats

    Example:
      >>> to_floats([1, 2, 3])
      [1.0, 2.0, 3.0]
    """

    return [float(item) for item in iterable_obj]


def arrays_are_close(a_array, b_array, rtol=1e-012, atol=1e-12, equal_nan=False, equal_inf=False):
    """
    Check for equivalence between two numpy arrays.

    :param a_array: numpy array
    :param b_array: numpy array
    :param rtol: relative tolerance
    :param atol: absolute tolerance
    :param equal_nan: consider nan values equivalent or not
    :param equal_inf: consider inf values equivalent or not
    :return: Boolean

    Example:
      >>> arrays_are_close(array([1,2,3]), array([1,2,3]))
      True
      >>> arrays_are_close(array([[1,2,3], [4, 5, 6]]), array([1,2,3]))
      False
      >>> arrays_are_close(array([[1,2,3], [4,5,6]]), array([[1,2,3], [4,5,6]]))
      True
      >>> arrays_are_close(array([[1,2,np.nan], [4,5,6]]), array([[1,2,np.nan], [4,5,6]]))
      False
      >>> arrays_are_close(array([[1,2,np.nan], [4,5,6]]), array([[1,2,np.nan], [4,5,6]]), equal_nan=True)
      True
    """
    if a_array.shape != b_array.shape:
        return False

    are_equal = []
    for a, b in np.nditer([a_array, b_array]):
        are_equal.append(areClose(a.item(0), b.item(0), rtol=rtol, atol=atol, equal_nan=equal_nan, equal_inf=equal_inf))

    return all(are_equal)


def point_solution(a_array, b_array):
    """
    finds a non-unique solution
    for a set of linear equations
    """

    try:
        return np.linalg.lstsq(a_array, b_array, rcond=None)[0]
    except:
        return None, None, None


def xyz_svd(xyz_array):
    """
    Calculates the SVD solution given a Numpy array.

    # modified after: 
    # http://stackoverflow.com/questions/15959411/best-fit-plane-algorithms-why-different-results-solved
    """

    try:
        result = np.linalg.svd(xyz_array)
    except:
        result = None

    return dict(result=result)


if __name__ == "__main__":

    import doctest
    doctest.testmod()
