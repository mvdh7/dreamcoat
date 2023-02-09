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
    route_line[-1] = stop
    for i in range(1, num - 1):
        route_line[i] = gcc.intermediate_point(start, stop, fraction=i / (num - 1))
    return route_line.T


def linspace_gc_route(route, num_per_segment=50, endpoint=True):
    """Apply `linspace_gc()` across a series of route points.

    Parameters
    ----------
    route : (array_like, array_like)
        The (longitude, latitude) of the route waypoints.
    num_per_segment : int, optional
        Number of points between each route point, by default 50
    endpoint : bool, optional
        Whether to include the final endpoint, by default True

    Returns
    -------
    (array_like, array_like)
        The (longitude, latitude) of the interpolated route.
    """
    route_line = []
    for i in range(1, route.shape[0]):
        route_line.append(
            linspace_gc(
                route[i - 1],
                route[i],
                num=num_per_segment,
                endpoint=i == route.shape[0] - 1 and endpoint,
            )
        )
    route_line = np.concatenate(route_line, axis=1)
    return route_line
