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


def cartesian_to_polar(x, y):
    """Convert from Cartesian to polar co-ordinates.

    Parameters
    ----------
    x : float
        x-axis value in Cartesian co-ordinates.
    y : float
        y-axis value in Cartesian co-ordinates.

    Returns
    -------
    theta
        Angle in polar co-ordinates.
    rho
        Radius in polar co-ordinates.
    """
    rho = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return theta, rho


def polar_to_cartesian(theta, rho):
    """Convert from polar to Cartesian co-ordinates.

    Parameters
    ----------
    theta
        Angle in polar co-ordinates.
    rho
        Radius in polar co-ordinates.

    Returns
    -------
    x : float
        x-axis value in Cartesian co-ordinates.
    y : float
        y-axis value in Cartesian co-ordinates.
    """
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return x, y
