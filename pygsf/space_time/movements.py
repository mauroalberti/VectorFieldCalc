# -*- coding: utf-8 -*-


from ..spatial.vectorial.vectorial import Point
from ..spatial.rasters.geoarray import GeoArray


def point_velocity(geoarray: GeoArray, currpoint: Point):
    """
    return the velocity components of a point in a velocity field

    Examples:
    """

    x, y, _ = currpoint.toXYZ()
    vx = geoarray.interpolate_bilinear(x, y, 0)
    vy = geoarray.interpolate_bilinear(x, y, 1)

    return vx, vy


def interpolate_RKF(geoarray: GeoArray, delta_time: Number, curr_Pt: Point):
    """
    Interpolate points according to RKF method.

    Examples:
    """

    K1_vx, K1_vy = point_velocity(geoarray,  curr_Pt)
    if K1_vx is None or K1_vy is None:
        return None, None

    K2_Pt = Point(
        curr_Pt.x + (0.25) * delta_time * K1_vx,
        curr_Pt.y + (0.25) * delta_time * K1_vy)

    K2_vx, K2_vy = point_velocity(geoarray,  K2_Pt)
    if K2_vx is None or K2_vy is None:
        return None, None

    K3_Pt = Point(
        curr_Pt.x + (3.0 / 32.0) * delta_time * K1_vx + (9.0 / 32.0) * delta_time * K2_vx,
        curr_Pt.y + (3.0 / 32.0) * delta_time * K1_vy + (9.0 / 32.0) * delta_time * K2_vy)

    K3_vx, K3_vy = point_velocity(geoarray,  K3_Pt)
    if K3_vx is None or K3_vy is None:
        return None, None

    K4_Pt = Point(
        curr_Pt.x + (1932.0 / 2197.0) * delta_time * K1_vx - (7200.0 / 2197.0) * delta_time * K2_vx + (7296.0 / 2197.0) * delta_time * K3_vx,
        curr_Pt.y + (1932.0 / 2197.0) * delta_time * K1_vy - (7200.0 / 2197.0) * delta_time * K2_vy + (7296.0 / 2197.0) * delta_time * K3_vy)

    K4_vx, K4_vy = point_velocity(geoarray,  K4_Pt)
    if K4_vx is None or K4_vy is None:
        return None, None

    K5_Pt = Point(
        curr_Pt.x + (439.0 / 216.0) * delta_time * K1_vx - (8.0) * delta_time * K2_vx + (3680.0 / 513.0) * delta_time * K3_vx - (845.0 / 4104.0) * delta_time * K4_vx,
        curr_Pt.y + (439.0 / 216.0) * delta_time * K1_vy - (8.0) * delta_time * K2_vy + (3680.0 / 513.0) * delta_time * K3_vy - (845.0 / 4104.0) * delta_time * K4_vy)

    K5_vx, K5_vy = point_velocity(geoarray,  K5_Pt)
    if K5_vx is None or K5_vy is None:
        return None, None

    K6_Pt = Point(
        curr_Pt.x - (8.0 / 27.0) * delta_time * K1_vx + (2.0) * delta_time * K2_vx - (3544.0 / 2565.0) * delta_time * K3_vx + (1859.0 / 4104.0) * delta_time * K4_vx - (
                          11.0 / 40.0) * delta_time * K5_vx,
        curr_Pt.y - (8.0 / 27.0) * delta_time * K1_vy + (2.0) * delta_time * K2_vy - (3544.0 / 2565.0) * delta_time * K3_vy + (1859.0 / 4104.0) * delta_time * K4_vy - (
                          11.0 / 40.0) * delta_time * K5_vy)

    K6_vx, K6_vy = point_velocity(geoarray,  K6_Pt)
    if K6_vx is None or K6_vy is None:
        return None, None

    rkf_4o_x = curr_Pt.x + delta_time * (
            (25.0 / 216.0) * K1_vx + (1408.0 / 2565.0) * K3_vx + (2197.0 / 4104.0) * K4_vx - (
            1.0 / 5.0) * K5_vx)
    rkf_4o_y = curr_Pt.y + delta_time * (
            (25.0 / 216.0) * K1_vy + (1408.0 / 2565.0) * K3_vy + (2197.0 / 4104.0) * K4_vy - (
            1.0 / 5.0) * K5_vy)
    temp_Pt = Point(
        rkf_4o_x,
        rkf_4o_y)

    interp_Pt_x = curr_Pt.x + delta_time * (
            (16.0 / 135.0) * K1_vx + (6656.0 / 12825.0) * K3_vx + (28561.0 / 56430.0) * K4_vx - (
            9.0 / 50.0) * K5_vx + (2.0 / 55.0) * K6_vx)
    interp_Pt_y = curr_Pt.y + delta_time * (
            (16.0 / 135.0) * K1_vy + (6656.0 / 12825.0) * K3_vy + (28561.0 / 56430.0) * K4_vy - (
            9.0 / 50.0) * K5_vy + (2.0 / 55.0) * K6_vy)
    interp_Pt = Point(
        interp_Pt_x,
        interp_Pt_y)

    interp_PT_error_estimate = interp_Pt.distance(temp_Pt)

    return interp_Pt, interp_PT_error_estimate

