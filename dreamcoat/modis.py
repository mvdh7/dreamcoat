import warnings
from datetime import date, datetime, timedelta
import xarray as xr, numpy as np

satellites = {"AQUA": "A", "TERRA": "T"}
variables = {"pic": "PIC", "poc": "POC"}


def _get_url_daily(satellite, year, month, day, resolution="9km"):
    return (
        "http://oceandata.sci.gsfc.nasa.gov/opendap/MODIS{}/L3SMI/".format(
            satellites[satellite.upper()]
        )
        + (
            "{year:4.0f}/{month:02.0f}{day:02.0f}/"
            + "{satellite}_MODIS.{year:4.0f}{month:02.0f}{day:02.0f}"
        ).format(year=year, month=month, day=day, satellite=satellite.upper())
        + ".L3m.DAY.{var_big}.{var_small}.{resolution}.NRT.nc".format(
            var_big="PIC", var_small="pic", resolution=resolution
        )
    )


def _get_url_8day(satellite, year, month, day, month_end, day_end):
    return (
        "http://oceandata.sci.gsfc.nasa.gov/opendap/MODIS{}/L3SMI/".format(
            satellites[satellite.upper()]
        )
        + (
            "{year:4.0f}/{month:02.0f}{day:02.0f}/"
            + "{satellite}_MODIS.{year:4.0f}{month:02.0f}{day:02.0f}_"
            + "{year:4.0f}{month_end:02.0f}{day_end:02.0f}"
        ).format(
            year=year,
            month=month,
            day=day,
            month_end=month_end,
            day_end=day_end,
            satellite=satellite.upper(),
        )
        + ".L3m.8D.{var_big}.{var_small}.4km.NRT.nc".format(
            var_big="PIC", var_small="pic"
        )
    )


def get_single_day(
    year,
    month,
    day,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    resolution="9km",
):
    """Download a single day of MODIS PIC data from both the Aqua and Terra satellites,
    combining the results into a single xarray Dataset.

    Parameters
    ----------
    year : int
        The year to download data for.
    month : int
        The month to download data for.
    day : int
        The day to download data for.
    latitude_min : float, optional
        Minimum latitude to subset data, by default -90.
    latitude_max : float, optional
        Maximum latitude to subset data, by default 90.
    longitude_min : float, optional
        Minimum longitude to subset data, by default -180.
    longitude_max : float, optional
        Maximum longitude to subset data, by default 180.
    resolution : str, optional
        Resolution to download, can be "4km" or "9km", by default "9km".

    Returns
    -------
    xarray.Dataset
        The combined PIC dataset with units mol/m**3.
    """
    modis_a = xr.open_dataset(
        _get_url_daily("aqua", year, month, day, resolution=resolution)
    )
    modis_t = xr.open_dataset(
        _get_url_daily("terra", year, month, day, resolution=resolution)
    )
    ix_lat = modis_a.lat.data[
        (modis_a.lat.data >= latitude_min) & (modis_a.lat.data <= latitude_max)
    ]
    ix_lon = modis_a.lon.data[
        (modis_a.lon.data >= longitude_min) & (modis_a.lon.data <= longitude_max)
    ]
    modis_a = modis_a.sel(lat=ix_lat, lon=ix_lon)
    modis_t = modis_t.sel(lat=ix_lat, lon=ix_lon)
    modis = modis_a.copy()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        modis["pic"] = (
            ("lat", "lon"),
            np.nanmean(np.array([modis_a.pic.data, modis_t.pic.data]), axis=0),
        )
    return modis


def get_days(
    date_min=None,
    date_max=None,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    resolution="9km",
):
    """Download a series of days of MODIS PIC data, combining both the Aqua and Terra
    satellites, and compile them into a single xarray Dataset.

    Parameters
    ----------
    date_min : str, optional
        First date of data to download in '%Y-%m-%d' format, by default None, in which
        case yesterday is used.
    date_max : str, optional
        Last date of data to download in '%Y-%m-%d' format, by default None, in which
        case yesterday is used.
    latitude_min : float, optional
        Minimum latitude to subset data, by default -90.
    latitude_max : float, optional
        Maximum latitude to subset data, by default 90.
    longitude_min : float, optional
        Minimum longitude to subset data, by default -180.
    longitude_max : float, optional
        Maximum longitude to subset data, by default 180.
    resolution : str, optional
        Resolution to download, can be "4km" or "9km", by default "9km".

    Returns
    -------
    xarray.Dataset
        The combined PIC dataset with units mol/m**3.
    """
    if date_min is None:
        date_min = (date.today() - timedelta(days=1)).strftime("%Y-%m-%d")
    if date_max is None:
        date_max = (date.today() - timedelta(days=1)).strftime("%Y-%m-%d")
    date_start = datetime(*[int(x) for x in date_min.split("-")])
    date_end = datetime(*[int(x) for x in date_max.split("-")])
    days_before = (date_end - date_start).days
    modis = []
    for d in range(days_before + 1):
        date_now = date_start + timedelta(days=d)
        print("Getting MODIS data for {}...".format(date_now.strftime("%Y-%m-%d")))
        modis.append(
            get_single_day(
                date_now.year,
                date_now.month,
                date_now.day,
                latitude_min=latitude_min,
                latitude_max=latitude_max,
                longitude_min=longitude_min,
                longitude_max=longitude_max,
                resolution=resolution,
            )
        )
        modis[d] = modis[d].expand_dims(dim={"date": [np.datetime64(date_now)]})
    modis = xr.concat(modis, "date")
    print("Got all MODIS data!")
    return modis


def _get_8day_end(date_start):
    date_plus7 = date_start + timedelta(days=7)
    if date_plus7.year == date_start.year:
        month_end = date_plus7.month
        day_end = date_plus7.day
    else:
        month_end = 12
        day_end = 31
    return month_end, day_end


def get_8day(
    year,
    month,
    day,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
):
    date_start = datetime(year, month, day)
    d = 0
    while d < 8:
        date_here = date_start - timedelta(days=d)
        month_end, day_end = _get_8day_end(date_here)
        try:
            modis_a = xr.open_dataset(
                get_modis_url_8day(
                    "aqua",
                    date_here.year,
                    date_here.month,
                    date_here.day,
                    month_end,
                    day_end,
                )
            )
            print(
                "8-day file found for {}-{:02.0f}-{:02.0f}!".format(
                    date_here.year, date_here.month, date_here.day
                )
            )
            d = 10
            break
        except OSError:
            print(
                (
                    "No 8-day file found for {}-{:02.0f}-{:02.0f}, "
                    + "trying one day earlier..."
                ).format(date_here.year, date_here.month, date_here.day)
            )
            d += 1
    try:
        modis_t = xr.open_dataset(
            get_modis_url_8day(
                "terra",
                date_here.year,
                date_here.month,
                date_here.day,
                month_end,
                day_end,
            )
        )
    except OSError:
        modis_t = modis_a.copy()
    ix_lat = modis_a.lat.data[
        (modis_a.lat.data >= latitude_min) & (modis_a.lat.data <= latitude_max)
    ]
    ix_lon = modis_a.lon.data[
        (modis_a.lon.data >= longitude_min) & (modis_a.lon.data <= longitude_max)
    ]
    modis_a = modis_a.sel(lat=ix_lat, lon=ix_lon)
    modis_t = modis_t.sel(lat=ix_lat, lon=ix_lon)
    modis = modis_a.copy()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        modis["pic"] = (
            ("lat", "lon"),
            np.nanmean(np.array([modis_a.pic.data, modis_t.pic.data]), axis=0),
        )
    return modis
