"""
dreamcoat.maps
==============
Geographic mapping and route planning.

Functions
---------
build_vptree
    Build vantage-point tree for Natural Earth Data coastline.
get_distance
    Calculate distance between a pair of points.
get_distance_gcc
    Calculate distance between a pair of points in km using great_circle_calculator.
get_distance_geodesic
    Calculate geodesic distance between a single pair of points in km.
get_route_distance
    Calculate distance from point to point along a route.
linspace_gc
    Return evenly spaced points along a great circle.
linspace_gc_waypoints
    Apply `linspace_gc` across a series of waypoints, aiming for evenly spaced points.
linspace_gc_waypoints_cs
    Apply `linspace_gc` across a series of waypoints, with a constant number of points
    between each pair.
map_point_to_route
    Find the closest point on a route to a new point.
map_point_to_section
    Find the closest point on a single line section to a given point.
"""

import itertools
from collections import namedtuple
import numpy as np
from geographiclib.geodesic import Geodesic
import great_circle_calculator.great_circle_calculator as gcc
from vptree import VPTree
import cartopy.io.shapereader as shpreader
from shapely import geometry
from . import degrees_decimal_minutes


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
    wp_interp = np.full((num - (not endpoint), 2), np.nan)
    wp_interp[0] = start
    if endpoint:
        wp_interp[-1] = stop
    for i in range(1, num - 1):
        wp_interp[i] = gcc.intermediate_point(start, stop, fraction=i / (num - 1))
    return wp_interp.T


def linspace_gc_waypoints(waypoints, num_approx=100, endpoint=True):
    """Apply `linspace_gc()` across a series of waypoints, aiming for approximately
    `num_approx` roughly evenly spaced final points.

    Parameters
    ----------
    waypoints : (array_like, array_like)
        The (longitude, latitude) of the waypoints.
    num_approx : int, optional
        Approximate total number of points in the final interpolation, by default 100.
    endpoint : bool, optional
        Whether to include the final endpoint, by default True.

    Returns
    -------
    (array_like, array_like)
        The (longitude, latitude) of the interpolated waypoints.
    """
    # Calculate relative distance for each route step
    r_distance = get_route_distance(waypoints)
    r_fraction = np.diff(r_distance) / r_distance[-1]
    # Divide num_approx into the steps
    r_pts = (num_approx - 1) * r_fraction
    r_pts_round = np.round(r_pts).astype(int)
    r_zero = r_pts_round == 0
    r_pts_round[r_zero] += 1
    total_points = np.sum(r_pts_round) + 1
    # Do individual section interpolations and concatenate them
    wp_interp = []
    for i, npts in enumerate(r_pts_round):
        wp_interp.append(
            linspace_gc(
                waypoints[:, i],
                waypoints[:, i + 1],
                num=npts + 1,
                endpoint=i == len(r_pts_round) - 1 and endpoint,
            )
        )
    wp_interp = np.concatenate(wp_interp, axis=1)
    return wp_interp


def linspace_gc_waypoints_cs(waypoints, num_per_segment=50, endpoint=True):
    """Apply `linspace_gc()` across a series of waypoints with a constant number
    of interpolated points between each pair of waypoints.

    Parameters
    ----------
    waypoints : (array_like, array_like)
        The (longitude, latitude) of the waypoints.
    num_per_segment : int, optional
        Number of points between each pair of waypoints, by default 50.
    endpoint : bool, optional
        Whether to include the final endpoint, by default True.

    Returns
    -------
    (array_like, array_like)
        The (longitude, latitude) of the interpolated route.
    """
    wp_interp = []
    for i in range(1, route.shape[0]):
        wp_interp.append(
            linspace_gc(
                route[:, i - 1],
                route[:, i],
                num=num_per_segment,
                endpoint=i == route.shape[0] - 1 and endpoint,
            )
        )
    wp_interp = np.concatenate(wp_interp, axis=1)
    return wp_interp


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
    """Calculate distance between a pair of points in km using either
    great_circle_calculator (default) or geographiclib (more accurate but ~50 times
    slower).

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
    if isinstance(lon_lat_1, degrees_decimal_minutes.LatLon):
        lon_lat_1 = (lon_lat_1.longitude_dd, lon_lat_1.latitude_dd)
    if isinstance(lon_lat_2, degrees_decimal_minutes.LatLon):
        lon_lat_2 = (lon_lat_2.longitude_dd, lon_lat_2.latitude_dd)
    return _get_distance_func(method)(lon_lat_1, lon_lat_2)


def get_heading(lon_lat_1, lon_lat_2):
    """Calculate heading from lon_lat_1 to lon_lat_2.

    Parameters
    ----------
    lon_lat_1 : array_like
        (longitude, latitude) of the first point in decimal degrees.
    lon_lat_2 : array_like
        (longitude, latitude) of the second point in decimal degrees.

    Returns
    -------
    float
        Heading in ° (with N as 0, E as 90, ...).
    """
    if isinstance(lon_lat_1, degrees_decimal_minutes.LatLon):
        lon_lat_1 = (lon_lat_1.longitude_dd, lon_lat_1.latitude_dd)
    if isinstance(lon_lat_2, degrees_decimal_minutes.LatLon):
        lon_lat_2 = (lon_lat_2.longitude_dd, lon_lat_2.latitude_dd)
    heading = gcc.bearing_at_p1(lon_lat_1, lon_lat_2)
    if heading < 0:
        heading += 360
    return heading


def get_route_distance(waypoints, method="gcc"):
    """Calculate distance from point to point along a route.

    Parameters
    ----------
    waypoints : array_like
        The route (longitude, latitude) waypoints as a size (2, n) np.array.
    method : str, optional
        Which method to use: 'gcc' (default) or 'geodesic'.

    Returns
    -------
    array_like
        Distance from the start of the route in km.
    """
    route_distance = np.full_like(waypoints[0], 0)
    for i in range(1, len(route_distance)):
        route_distance[i] = (
            get_distance(waypoints[:, i - 1], waypoints[:, i], method=method)
            + route_distance[i - 1]
        )
    return route_distance


SectionPoint = namedtuple(
    "SectionPoint",
    ("distance_on_section", "distance_perpendicular", "lon_lat", "in_section"),
)
RoutePoint = namedtuple(
    "RoutePoint",
    ("distance_on_route", "lon_lat", "in_route", "waypoint_nearest"),
)


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
    return SectionPoint(
        section_distance_pt,
        distance_perpendicular,
        lon_lat_pt_section,
        in_section,
    )


def map_point_to_route(route, lon_lat, extrapolate=False, verbose=False):
    """Find the closest point on a `route` to a new `lon_lat` point.

    Parameters
    ----------
    route : Route
        A Route object containing the route's waypoints.
    lon_lat : array_like
        (longitude, latitude) of the new off-route point.
    extrapolate : bool, optional
        Whether to extrapolate beyond the ends of the route, by default False.
    verbose : bool, optional
        Whether to report on progress, by default False.

    Results
    -------
    namedtuple (RoutePoint) containing fields
        distance_on_route : float
            Distance along the route up to the mapping of `lon_lat_pt` onto the route
            in km.
        lon_lat : (float, float)
            (longitude, latitude) of the mapping of `lon_lat` onto the route in decimal
            degrees.
        in_route : bool
            Whether the mapping of `lon_lat` falls within the route.
        waypoint_nearest : int
            The index of the nearest waypoint in `route.waypoints`.
    """

    def printv(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    nearest = route.find_nearest_waypoint(lon_lat)
    section = None
    # Get the next section after the nearest waypoint, if possible
    section_next = None
    if nearest < route.size - 1:
        printv("Nearest waypoint is not the last waypoint.")
        section_next = map_point_to_section(
            route.waypoints[:, nearest],
            route.waypoints[:, nearest + 1],
            lon_lat,
            method=route.distance_method,
        )
    # Get the section directly before the nearest waypoint, if possible
    section_prev = None
    if nearest > 0:
        printv("Nearest waypoint is not the first waypoint.")
        section_prev = map_point_to_section(
            route.waypoints[:, nearest],
            route.waypoints[:, nearest - 1],
            lon_lat,
            method=route.distance_method,
        )
        # If using previous section, distance needs to be negative
        section_prev = SectionPoint(
            -section_prev.distance_on_section,
            section_prev.distance_perpendicular,
            section_prev.lon_lat,
            section_prev.in_section,
        )
    # Decide which of the two sections to use
    if section_next is not None and section_prev is not None:
        printv("Nearest waypoint is neither the first nor the last waypoint.")
        if section_next.in_section and section_prev.in_section:
            printv(
                " Point falls within both sections either side of the nearest waypoint."
            )
            if (
                section_next.distance_perpendicular
                <= section_prev.distance_perpendicular
            ):
                printv("  Point is closer to the section after the nearest waypoint.")
                section = section_next
            else:
                printv("  Point is closer to the section before the nearest waypoint.")
                section = section_prev
        elif section_next.in_section:
            printv(" Point falls only in the section after the nearest waypoint.")
            section = section_next
        elif section_prev.in_section:
            printv(" Point falls only in the section before the nearest waypoint.")
            section = section_prev
        else:
            printv(" Point does not fall in either section.")
            section = SectionPoint(0, 0, route.waypoints[:, nearest], True)
    elif section_next is not None:
        printv("Nearest waypoint is the first waypoint.")
        if section_next.in_section:
            printv(" Point falls within the first section.")
            section = section_next
        else:
            printv(" Point does not fall within the first section.")
            if extrapolate:
                printv("  Extrapolating before the first section.")
                section = section_next
            else:
                printv("  Returning the first waypoint.")
                section = SectionPoint(0, 0, route.waypoints[:, nearest], True)
    elif section_prev is not None:
        printv("Nearest waypoint is the last waypoint.")
        if section_prev.in_section:
            printv(" Point falls within the last section.")
            section = section_prev
        else:
            printv(" Point does not fall within the last section.")
            if extrapolate:
                printv("  Extrapolating after the last section.")
                section = section_prev
            else:
                printv("  Returning the last waypoint.")
                section = SectionPoint(0, 0, route.waypoints[:, nearest], True)
    # Convert section properties to route properties
    route_point = RoutePoint(
        route.distance[nearest] + section.distance_on_section,
        section.lon_lat,
        section.in_section,
        nearest,
    )
    return route_point


class Route:
    """A set of waypoints defining a route, which is assumed to travel in straight lines
    (i.e., along great circles) between the waypoints.

    Parameters
    ----------
    waypoints : np.ndarray
        A (2, n) NumPy array of the (longitude, latitude) values of the waypoints.
    distance_method : str, optional
        Which method is used to calculate distances, either
          "gcc" (default) - uses great_circle_calculator
          "geodesic" - uses geographiclib.geodesic.Geodesic.WGS84
        The first is an order of magnitude faster but a bit less accurate.

    Attributes
    ----------
    waypoints : np.ndarray
        The (2, n) NumPy array of the (longitude, latitude) values of the waypoints.
    size : int
        The number of waypoints (i.e., n).
    distance_method : str
        The method being used to calculate distances ("gcc" or "geodesic").
    distance : np.ndarray
        The distance of each waypoint along the route in km.
    vpt : vptree.VPTree
        A vantage-point tree for the route.
    wp_interp : np.ndarray
        A (2, n) NumPy array of the (longitude, latitude) values of a higher-resolution
        set of points interpolated between the waypoints.  Appears after running the
        interp method.
    wp_distance : np.ndarray
        The distance of each interpolated point in wp_interp along the route in km.
        Appears after running the interp method.

    Methods
    -------
    build_vpt
        Create a vantage-point tree for the route.
    find_nearest_waypoint
        Find the index of the waypoint closest to a given new point.
    get_distance
        Find the distance along the route corresponding to the closest point on the
        route to a new point.
    get_lon_lat
        Calculate the longitude and latitude of a given distance along the route.
    interp
        Interpolate the route to a higher resolution.
    map_point_to_route
        Find the closest point on the route to a new point.
    """

    def __init__(self, waypoints, distance_method="gcc"):
        assert (
            isinstance(waypoints, np.ndarray)
            and waypoints.shape[0] == 2
            and waypoints.ndim == 2
        ), "`waypoints` must be a numpy array of shape (2, n)."
        assert np.all(
            (waypoints[0] >= -180) & (waypoints[0] <= 180)
        ), "All longitudes must be from -180 to 180."
        assert np.all(
            (waypoints[1] >= -90) & (waypoints[1] <= 90)
        ), "All latitudes must be from -90 to 90."
        self.waypoints = waypoints
        self.size = waypoints.shape[1]
        self.distance_method = distance_method
        self.distance = get_route_distance(waypoints, method=distance_method)
        self.build_vpt()

    def build_vpt(self):
        """Create a vantage-point tree that can be used to find nearest neighbours."""
        self.vpt = VPTree(self.waypoints.T, _get_distance_func(self.distance_method))

    def find_nearest_waypoint(self, lon_lat):
        """Find the index of the waypoint closest to a new point.

        Parameters
        ----------
        lon_lat : array-like
            A two-element list or array containing the (longitude, latitude) of a new
            point, not necessarily on the route.

        Returns
        -------
        int
            The index of the waypoint closest to the new point.
        """
        nearest = self.vpt.get_nearest_neighbor(lon_lat)[1]
        return (
            (self.waypoints[0] == nearest[0]) & (self.waypoints[1] == nearest[1])
        ).argmax()

    def map_point_to_route(self, lon_lat, extrapolate=False, verbose=False):
        """Find the closest point on the route to a new point.

        Parameters
        ----------
        lon_lat : array-like
            (longitude, latitude) of the new off-route point.
        extrapolate : bool, optional
            Whether to extrapolate beyond the ends of the route, by default False.
        verbose : bool, optional
            Whether to report on progress, by default False.

        Results
        -------
        namedtuple (RoutePoint) containing fields
            distance_on_route : float
                Distance along the route up to the mapping of `lon_lat` onto the route
                in km.
            lon_lat : (float, float)
                (longitude, latitude) of the mapping of `lon_lat` onto the route in
                decimal degrees.
            in_route : bool
                Whether the mapping of `lon_lat` falls within the route.
            waypoint_nearest : int
                The index of the nearest waypoint.
        """
        return map_point_to_route(
            self, lon_lat, extrapolate=extrapolate, verbose=verbose
        )

    def interp(self, num_approx=100, endpoint=True):
        """Interpolate the route to a higher resolution, with the interpolated route
        stored in the attribute `wp_interp` and corresponding route distances in
        `wp_distance`.

        The interpolation runs such that the exact positions of the original waypoints
        are still included and the interpolated points are distributed as evenly as
        possible between the waypoints.

        Parameters
        ----------
        num_approx : int, optional
            Approximate total number of interpolated points, by default 100.
        endpoint : bool, optional
            Whether to include the very final waypoint in the interpolated set, by
            default True.
        """
        self.wp_interp = linspace_gc_waypoints(
            self.waypoints, num_approx=num_approx, endpoint=endpoint
        )
        self.wp_distance = get_route_distance(
            self.wp_interp, method=self.distance_method
        )

    def get_lon_lat(self, distance, extrapolate=False):
        """Get the longitude and latitude of a given distance along the route.

        Parameters
        ----------
        distance : float
            A distance along the route in km.
        extrapolate : bool, optional
            Whether to extrapolate beyond the ends of the route, if necessary; by
            default False.

        array-like
            (longitude, latitude) of the point along the route at the specified
            distance.
        """
        out_of_range = False
        if distance in self.distance:
            nearest = (self.distance == distance).argmax()
            lon_lat = self.waypoints[:, nearest]
        else:
            if distance < 0:
                after = 1
                before = 0
                out_of_range = True
            elif distance > self.distance.max():
                after = -1
                before = -2
                out_of_range = True
            else:
                after_point = distance < self.distance
                after = after_point.argmax()
                before = after - 1
            if not extrapolate and out_of_range:
                if distance < 0:
                    lon_lat = self.waypoints[:, 0]
                elif distance > self.distance.max():
                    lon_lat = self.waypoints[:, -1]
            else:
                section_fraction = (distance - self.distance[before]) / (
                    self.distance[after] - self.distance[before]
                )
                lon_lat = gcc.intermediate_point(
                    self.waypoints[:, before],
                    self.waypoints[:, after],
                    fraction=section_fraction,
                )
        return lon_lat

    def get_distance(self, lon_lat, extrapolate=False, verbose=False):
        """Get the route distance corresponding to a given longitude and latitude.

        Parameters
        ----------
        lon_lat : array-like
            (longitude, latitude) of the new off-route point.
        extrapolate : bool, optional
            Whether to extrapolate beyond the ends of the route, by default False.
        verbose : bool, optional
            Whether to report on progress, by default False.

        Returns
        -------
        float
            Distance along the route in km.
        """
        rp = self.map_point_to_route(lon_lat, extrapolate=extrapolate, verbose=verbose)
        return rp.distance_on_route


def geodesic_distance(lon_lat_1, lon_lat_2):
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
        / 1000
    )


def coastline_coords(lat_range, lon_range, resolution="10m"):
    """Get list of coastline coordinates from Natural Earth Data coastline sections
    that intersect with the lat_range and lon_range.

    Parameters
    ----------
    lat_range : array_like
        (min, max) latitude in decimal degrees.
    lon_range : array_like
        (min, max) longitude in decimal degrees.
    resolution : str, optional
        Natural Earth Data resolution, by default "10m" (can also be "50m" or "110m";
        default "10m" is the most detailed).

    Returns
    -------
    list
        Coordinates of the matching coastlines.
    """
    assert resolution in ["10m", "50m", "110m"]
    x = shpreader.natural_earth(
        resolution=resolution, category="physical", name="coastline"
    )
    cl = shpreader.Reader(x)
    coastlines = cl.records()
    my_bounds = geometry.Polygon(
        [
            (lon_range[0], lat_range[0]),
            (lon_range[0], lat_range[1]),
            (lon_range[1], lat_range[1]),
            (lon_range[1], lat_range[0]),
        ]
    )
    in_bounds = [cl for cl in coastlines if cl.geometry.intersects(my_bounds)]
    coords = list(itertools.chain(*[ib.geometry.coords for ib in in_bounds]))
    return coords


def build_vptree(lat_range, lon_range, **kwargs):
    """Build vantage-point tree from Natural Earth Data coastline sections that
    intersect with the lat_range and lon_range.

    Once constructed, use `vpt.get_nearest_neighbor((lon, lat))` to find the geodesic
    distance of nearest point on the coastline to (lon, lat) in km.

    Based on https://stackoverflow.com/a/45480555

    Parameters
    ----------
    lat_range : array_like
        (min, max) latitude in decimal degrees.
    lon_range : array_like
        (min, max) longitude in decimal degrees.

    Returns
    -------
    vpt : vptree.VPTree
        Vantage-point tree containing the matching coastline sections.
    """
    coords = coastline_coords(lat_range, lon_range, **kwargs)
    return VPTree(coords, geodesic_distance)
