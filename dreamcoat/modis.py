from datetime import datetime, timedelta
import xarray as xr, numpy as np

satellites = {"AQUA": "A", "TERRA": "T"}


def _get_modis_url_daily(satellite, year, month, day, resolution="9km"):
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


def _get_modis_url_8day(satellite, year, month, day, month_end, day_end):
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


def get_modis_daily(
    year,
    month,
    day,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    resolution="9km",
):
    modis_a = xr.open_dataset(
        _get_modis_url_daily("aqua", year, month, day, resolution=resolution)
    )
    modis_t = xr.open_dataset(
        _get_modis_url_daily("terra", year, month, day, resolution=resolution)
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
    modis["pic"] = (
        ("lat", "lon"),
        np.nanmean(np.array([modis_a.pic.data, modis_t.pic.data]), axis=0),
    )
    return modis


def get_modis_days_before(
    year,
    month,
    day,
    days_before=7,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    resolution="9km",
):
    date_end = datetime(year, month, day)
    date_start = date_end - timedelta(days=days_before)
    modis = []
    for d in range(days_before + 1):
        date = date_start + timedelta(days=d)
        print("Getting MODIS data for {}...".format(date.strftime("%Y-%m-%d")))
        modis.append(
            get_modis_daily(
                date.year,
                date.month,
                date.day,
                latitude_min=latitude_min,
                latitude_max=latitude_max,
                longitude_min=longitude_min,
                longitude_max=longitude_max,
                resolution=resolution,
            )
        )
        modis[d] = modis[d].expand_dims(dim={"date": [np.datetime64(date)]})
    modis = xr.concat(modis, "date")
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


def get_modis_8day(
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
    modis["pic"] = (
        ("lat", "lon"),
        np.nanmean(np.array([modis_a.pic.data, modis_t.pic.data]), axis=0),
    )
    return modis
