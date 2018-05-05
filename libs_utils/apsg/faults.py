# -*- coding: utf-8 -*-


from .exceptions import *


def rake_to_apsg_movsense(rake):
    """
    Convert rake to movement sense according to APSG convention.

    :param rake:
    :return:
    """

    if 0 < rake < 180:  # reverse faults according to Aki & Richards, 1980 convention
        return 1
    elif -180 < rake < 0:  # normal faults according to Aki & Richards, 1980 convention
        return -1
    elif abs(rake) == 0.0 or abs(rake) == 180.0:
        raise RakeInputException("Currently trascurrent data (rake = +/-180 or = 0.0) are not handled in plots")
    else:
        raise RakeInputException("Input rake value not acceptable")


def movsense_to_apsg_movsense(str_val):
    """
    Convert text movement sense convention to APSG convention.

    :param str_val:
    :return:
    """

    if str_val == "R":
        return 1
    elif str_val == "N":
        return -1
    else:
        raise RakeInputException("Input rake value not acceptable")

