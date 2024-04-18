"""
dreamcoat.convert
=================
Conversions between different units and formats.

Functions
---------
extent_to_nsew
    Get "N", "S", "E" and "W" from map extent values depending on their sign.
knots_to_kph
    Knots into kilometers per hour.
kph_to_knots
    Kilometers per hour into knots.
nm_to_km
    Nautical miles into kilometers.
km_to_nm
    Kilometers into nautical miles.
cartesian_to_polar
    Cartesian into polar coordinates.
polar_to_cartesian
    Polar into Cartesian coordinates.
"""

import numpy as np


def extent_to_nsew(map_extent):
    """Get "N", "S", "E" and "W" from map extent values depending on their sign.

    Parameters
    ----------
    map_extent : list of float
        [west, east, north, south] map extents in decimal degrees.

    Returns
    -------
    list of str
        ["E"/"W", "E"/"W", "N"/"S", "N"/"S"] as appropriate for the given extents.
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
    """Convert from Cartesian to polar co-ordinates, where theta is zero in the N
    position and increases clockwise.

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
    theta = np.arctan2(y, -x) - np.pi / 2
    return theta, rho


def polar_to_cartesian(theta, rho):
    """Convert from polar to Cartesian co-ordinates, where theta is zero in the N
    position and increases clockwise.

    Parameters
    ----------
    theta : float
        Angle in polar co-ordinates.
    rho : float
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
