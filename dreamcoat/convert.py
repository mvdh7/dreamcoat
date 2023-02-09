import os, shutil
import numpy as np, xarray as xr


def nc_to_zarr_zip(filename_nc, make_zip=True, delete_nc=False):
    """Convert netcdf to zarr and a zip archive of the zarr.

    Parameters
    ----------
    filename_nc : str
        File path and name for the netcdf file.
    make_zip : bool, optional
        Whether to make the zip archive, by default True.
    delete_nc : bool, optional
        Whether to delete the original netcdf file, by default False.
    """
    nc = xr.open_dataset(filename_nc)
    filename_zarr = ".zarr".join(filename_nc.rsplit(".nc", 1))
    nc.to_zarr(filename_zarr)
    if make_zip:
        shutil.make_archive(filename_zarr, "zip", filename_zarr)
    if delete_nc:
        os.remove(filename_nc)


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


def dd_to_ddm(dd):
    """Convert decimal degrees into degrees decimal minutes.

    Parameters
    ----------
    dd : float
        (longitude, latitude) in decimal degrees.

    Returns
    -------
    str
        The position(s) in degrees decimal minutes.
    """
    single = dd.ndim == 1
    if single:
        dd = [dd]
    ddm = []
    for d in dd:
        lon, lat = d
        ddm.append(
            "{:02.0f}°{:04.1f}'{}, {:03.0f}°{:04.1f}'{}".format(
                np.floor(np.abs(lat)),
                60 * (np.abs(lat) % 1),
                "SN"[int(lat > 0)],
                np.floor(np.abs(lon)),
                60 * (np.abs(lon) % 1),
                "WE"[int(lon > 0)],
            )
        )
    if single:
        ddm = ddm[0]
    return ddm
