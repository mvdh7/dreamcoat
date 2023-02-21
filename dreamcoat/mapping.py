from collections import namedtuple
import numpy as np
from geographiclib.geodesic import Geodesic
import great_circle_calculator.great_circle_calculator as gcc


def linspace_gc(start, stop, num=50, endpoint=True):
    """Return evenly spaced points along a great circle.

    Parameters
    ----------
    start : (float, float)
        The (longitude, latitude) of the start point.
    stop : (float, float)
        The (longitude, latitude) of the end point.
    num : int, optional
        Number of points including start and stop, by default 50.
    endpoint : bool, optional
        Whether to include the endpoint, by default True.

    Returns
    -------
    (array_like, array_like)
        The (longitude, latitude) of the interpolated route.
    """
    route_line = np.full((num - (not endpoint), 2), np.nan)
    route_line[0] = start
    if endpoint:
        route_line[-1] = stop
    for i in range(1, num - 1):
        route_line[i] = gcc.intermediate_point(start, stop, fraction=i / (num - 1))
    return route_line.T


def linspace_gc_route(route, num_approx=100, endpoint=True):
    """Apply `linspace_gc()` across a series of route points aiming for approximately
    `num_approx` roughly evenly spaced final points.

    Parameters
    ----------
    route : (array_like, array_like)
        The (longitude, latitude) of the route waypoints.
    num_approx : int, optional
        Approximate total number of points in the final route, by default 100.
    endpoint : bool, optional
        Whether to include the final endpoint, by default True.

    Returns
    -------
    (array_like, array_like)
        The (longitude, latitude) of the interpolated route.
    """
    # Calculate relative distance for each route step
    r_distance = get_route_distance(route)
    r_fraction = np.diff(r_distance) / r_distance[-1]
    # Divide num_approx into the steps
    r_pts = (num_approx - 1) * r_fraction
    r_pts_round = np.round(r_pts).astype(int)
    r_zero = r_pts_round == 0
    r_pts_round[r_zero] += 1
    total_points = np.sum(r_pts_round) + 1
    # Do individual section interpolations and concatenate them
    route_line = []
    for i, npts in enumerate(r_pts_round):
        route_line.append(
            linspace_gc(
                route[:, i],
                route[:, i + 1],
                num=npts + 1,
                endpoint=i == len(r_pts_round) - 1 and endpoint,
            )
        )
    route_line = np.concatenate(route_line, axis=1)
    return route_line


def linspace_gc_route_cs(route, num_per_segment=50, endpoint=True):
    """Apply `linspace_gc()` across a series of route points with a constant number
    of points between each route segment.

    Parameters
    ----------
    route : (array_like, array_like)
        The (longitude, latitude) of the route waypoints.
    num_per_segment : int, optional
        Number of points between each route point, by default 50.
    endpoint : bool, optional
        Whether to include the final endpoint, by default True.

    Returns
    -------
    (array_like, array_like)
        The (longitude, latitude) of the interpolated route.
    """
    route_line = []
    for i in range(1, route.shape[0]):
        route_line.append(
            linspace_gc(
                route[:, i - 1],
                route[:, i],
                num=num_per_segment,
                endpoint=i == route.shape[0] - 1 and endpoint,
            )
        )
    route_line = np.concatenate(route_line, axis=1)
    return route_line


def get_route_distance(route):
    """Calculate distance from point to point along a route.

    Parameters
    ----------
    route : (array_like, array_like)
        The route (longitude, latitude) points.

    Returns
    -------
    array_like
        Distance from the start of the route in km.
    """
    route_distance = np.full_like(route[0], 0)
    for i in range(1, len(route_distance)):
        route_distance[i] = (
            gcc.distance_between_points(route[:, i - 1], route[:, i])
            + route_distance[i - 1]
        )
    route_distance /= 1e3
    return route_distance


def get_distance_geodesic(lon_lat_1, lon_lat_2):
    """Calculate geodesic distance between a single pair of points in km, based on
    https://stackoverflow.com/a/45480555

    Parameters
    ----------
    lon_lat_1 : array_like
        (longitude, latitude) of the first point in decimal degrees.
    lon_lat_2 : array_like
        (longitude, latitude) of the second point in decimal degrees.

    Returns
    -------
    float
        Distance between the two points in km.
    """
    return (
        Geodesic.WGS84.Inverse(lon_lat_1[1], lon_lat_1[0], lon_lat_2[1], lon_lat_2[0])[
            "s12"
        ]
        / 1e3
    )


def get_distance_gcc(lon_lat_1, lon_lat_2):
    """Calculate distance between a pair of points in km using great_circle_calculator.

    Parameters
    ----------
    lon_lat_1 : array_like
        (longitude, latitude) of the first point in decimal degrees.
    lon_lat_2 : array_like
        (longitude, latitude) of the second point in decimal degrees.

    Returns
    -------
    float
        Distance between the two points in km.
    """
    return gcc.distance_between_points(lon_lat_1, lon_lat_2, unit="kilometers")


def _get_distance_func(method):
    assert method in ["gcc", "geodesic"], "`method` must be either 'gcc' or 'geodesic'"
    if method == "gcc":
        return get_distance_gcc
    elif method == "geodesic":
        return get_distance_geodesic


def get_distance(lon_lat_1, lon_lat_2, method="gcc"):
    """Calculate distance between a pair of points in km using either great_circle_calculator.
    (default) or geographiclib (more accurate but ~50 times slower).
    Parameters
    ----------
    lon_lat_1 : array_like
        (longitude, latitude) of the first point in decimal degrees.
    lon_lat_2 : array_like
        (longitude, latitude) of the second point in decimal degrees.
    method : str, optional
        Which method to use: 'gcc' (default) or 'geodesic'.

    Returns
    -------
    float
        Distance between the two points in km.
    """
    return _get_distance_func(method)(lon_lat_1, lon_lat_2)


def map_point_to_section(lon_lat_1, lon_lat_2, lon_lat_pt, method="gcc"):
    """Find the closest point on a line section defined by `lon_lat_1` and `lon_lat_2`
    to a given point `lon_lat_pt` and calculate various associated properties.

    Parameters
    ----------
    lon_lat_1 : array_like
        (longitude, latitude) of the first point on the section in decimal degrees.
    lon_lat_2 : array_like
        (longitude, latitude) of the second point on the section in decimal degrees.
    lon_lat_pt : array_like
        (longitude, latitude) of the off-section point in decimal degrees.
    method : str, optional
        Which method to use for calculating distances: 'gcc' (default) or 'geodesic'.

    Returns
    -------
    namedtuple (SectionPoint) containing fields
        distance_on_section : float
            Distance between `lon_lat_1` and the mapping of `lon_lat_pt` onto the
            section in km.
        distance_perpendicular : float
            Distance between `lon_lat_pt` and its mapping onto the section in km.
        lon_lat : (float, float)
            (longitude, latitude) of the mapping of `lon_lat_pt` onto the section in
            decimal degrees.
        in_section : bool
            Whether the mapping of `lon_lat_pt` falls within the section.
    """
    distance_func = _get_distance_func(method)
    bearing_section = np.deg2rad(gcc.bearing_at_p1(lon_lat_1, lon_lat_2))
    bearing_point = np.deg2rad(gcc.bearing_at_p1(lon_lat_1, lon_lat_pt))
    angle_to_section = bearing_point - bearing_section
    distance_r0_pt = distance_func(lon_lat_1, lon_lat_pt)
    section_distance = distance_func(lon_lat_1, lon_lat_2)
    section_distance_pt = distance_r0_pt * np.sin(np.deg2rad(90) - angle_to_section)
    lon_lat_pt_section = gcc.point_given_start_and_bearing(
        lon_lat_1, np.rad2deg(bearing_section), section_distance_pt, unit="kilometers"
    )
    in_section = section_distance_pt >= 0 and section_distance_pt <= section_distance
    distance_perpendicular = distance_func(lon_lat_pt, lon_lat_pt_section)
    point_on_section = namedtuple(
        "SectionPoint",
        ("distance_on_section", "distance_perpendicular", "lon_lat", "in_section"),
    )
    return point_on_section(
        section_distance_pt,
        distance_perpendicular,
        lon_lat_pt_section,
        in_section,
    )
