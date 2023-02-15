import numpy as np
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
