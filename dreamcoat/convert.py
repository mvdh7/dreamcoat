import numpy as np


def extent_to_nsew(map_extent):
    """Get N, S, E and W from map extent values depending on their sign.

    Parameters
    ----------
    map_extent : list
        [West, east, north, south] map extents in decimal degrees.

    Returns
    -------
    list
        [E/W, E/W, N/S, N/S] as appropriate for the given extents.
    """
    nsew = [
        " EW"[int(np.sign(map_extent[0]))],
        " EW"[int(np.sign(map_extent[1]))],
        " NS"[int(np.sign(map_extent[2]))],
        " NS"[int(np.sign(map_extent[3]))],
    ]
    return nsew


def knots_to_kph(knots):
    """Convert speed in knots to kilometers per hour.

    Parameters
    ----------
    knots : float
        Speed in knots.

    Returns
    -------
    float
        Speed in km/h.
    """
    return knots * 1.852


def kph_to_knots(kph):
    """Convert speed in kilometers per hour to knots.

    Parameters
    ----------
    kph : float
        Speed in km/h.

    Returns
    -------
    float
        Speed in knots.
    """
    return kph / 1.852


def nm_to_km(nm):
    """Convert distance in nautical miles to kilometers.

    Parameters
    ----------
    nm : float
        Distance in nautical miles.

    Returns
    -------
    float
        Distance in km.
    """
    return nm * 1.852


def km_to_nm(km):
    """Convert distance in kilometers to nautical miles.

    Parameters
    ----------
    km : float
        Distance in kilometers.

    Returns
    -------
    float
        Distance in nautical miles.
    """
    return km / 1.852
