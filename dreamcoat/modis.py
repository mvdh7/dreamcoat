import warnings, os
from datetime import date, datetime, timedelta
import xarray as xr, numpy as np
from . import meta

satellites = {"AQUA": "A", "TERRA": "T"}
variables = {"pic": "PIC", "poc": "POC"}


def _get_filename_daily(satellite, year, month, day, resolution="9km"):
    return (
        "{}_MODIS.".format(satellite.upper())
        + "{:4.0f}{:02.0f}{:02.0f}.L3m.DAY.".format(year, month, day)
        + "{var_big}.{var_small}.{resolution}.NRT".format(
            var_big="PIC", var_small="pic", resolution=resolution
        )
    )


def _get_url_daily(satellite, year, month, day, appkey=None, resolution="9km"):
    assert resolution in ["4km", "9km"], "Resolution must be either '4km' or '9km'."
    if appkey is None:
        appkey = meta.get_dat_data("nasa_appkey")
    return (
        "https://oceandata.sci.gsfc.nasa.gov/ob/getfile/"
        + _get_filename_daily(satellite, year, month, day, resolution=resolution)
        + ".nc?appkey={}".format(appkey)
    )


def download_single_day(
    satellite, year, month, day, appkey=None, filepath=".", resolution="9km"
):
    url = _get_url_daily(
        satellite, year, month, day, appkey=appkey, resolution=resolution
    )
    command = (
        "wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies"
        + "--auth-no-challenge=on --content-disposition {}".format(url)
    )
    if filepath is not None:
        command += " -P {}".format(filepath)
    os.system(command)


def get_single_day(
    year,
    month,
    day,
    appkey=None,
    delete_nc=False,
    filepath=".",
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
    appkey : str, optional
        The AppKey for your NASA Earthdata account, by default None, in which case it
        is taken from `.dreamcoat/nasa_appkey.dat`.
    delete_nc : bool, optional
        Whether to delete the nc files after converting, by default False.
    filepath : str, optional
        The file path where the data are to be saved, without a trailing separator,
        by default "".
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
    # Download the netcdf datasets from both satellites
    download_single_day(
        "aqua",
        year,
        month,
        day,
        appkey=appkey,
        filepath=filepath,
        resolution=resolution,
    )
    download_single_day(
        "terra",
        year,
        month,
        day,
        appkey=appkey,
        filepath=filepath,
        resolution=resolution,
    )
    # Import the dataset from each satellite separately
    filename_aqua = os.sep.join(
        (
            filepath,
            _get_filename_daily("aqua", year, month, day, resolution=resolution)
            + ".nc",
        )
    )
    filename_terra = os.sep.join(
        (
            filepath,
            _get_filename_daily("terra", year, month, day, resolution=resolution)
            + ".nc",
        )
    )
    modis_a = xr.open_dataset(filename_aqua)
    modis_t = xr.open_dataset(filename_terra)
    # Cut to the requested latitude and longitude ranges
    ix_lat = modis_a.lat.data[
        (modis_a.lat.data >= latitude_min) & (modis_a.lat.data <= latitude_max)
    ]
    ix_lon = modis_a.lon.data[
        (modis_a.lon.data >= longitude_min) & (modis_a.lon.data <= longitude_max)
    ]
    modis_a = modis_a.sel(lat=ix_lat, lon=ix_lon)
    modis_t = modis_t.sel(lat=ix_lat, lon=ix_lon)
    # Take the mean of the two satellites' data to get the final combined dataset
    modis = modis_a.copy()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        modis["pic"] = (
            ("lat", "lon"),
            np.nanmean(np.array([modis_a.pic.data, modis_t.pic.data]), axis=0),
        )
    if delete_nc:
        os.remove(filename_aqua)
        os.remove(filename_terra)
    return modis


def get_days(
    date_min=None,
    date_max=None,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    appkey=None,
    delete_nc=False,
    filepath=".",
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
    appkey : str, optional
        The AppKey for your NASA Earthdata account, by default None, in which case it
        is taken from `.dreamcoat/nasa_appkey.dat`.
    delete_nc : bool, optional
        Whether to delete the nc files after converting, by default False.
    filepath : str, optional
        The file path where the data are to be saved, without a trailing separator,
        by default "".
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
                appkey=appkey,
                delete_nc=delete_nc,
                filepath=filepath,
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


def open_days(
    filepath=".",
    date_start="1900-01-01",
    date_end="2300-01-01",
):
    """Open a set of downloaded MODIS zarr files in a single dataframe.

    Parameters
    ----------
    filepath : str, optional
        The file path where the MODIS zarr files can be found, by default ".".
    date_start : str, optional
        The earliest date to include data from in '%Y-%m-%d' format, by default
        "1900-01-01".
    date_end : str, optional
        The latest date to include data from in '%Y-%m-%d' format, by default
        "2300-01-01".
    """
    # Add separator to filepath if needed
    if not filepath.endswith(os.sep):
        filepath += os.sep
    # Get a list of all available MODIS files
    modis_files = [
        f
        for f in os.listdir(filepath)
        if f.startswith("modis_pic_") and f.endswith(".zarr")
    ]
    modis_files.sort()
    # Get the date associated with each file
    modis_dates = {
        f: np.datetime64(f.split("_")[-1].split(".")[0]) for f in modis_files
    }
    # Cut the list of files based on date_start and date_end
    modis_files = [
        f
        for f in modis_files
        if (modis_dates[f] >= np.datetime64(date_start))
        and (modis_dates[f] <= np.datetime64(date_end))
    ]
    # Import and concatenate all the files
    modis = xr.concat(
        [xr.open_dataset(filepath + f, engine="zarr") for f in modis_files], "date"
    )
    # Rename coordinates & dimensions for consistency with CMEMS
    modis = modis.rename({"lat": "latitude", "lon": "longitude", "date": "time"})
    modis = modis.sortby("latitude")
    # Convert units to mmol / m**3
    modis["pic"] *= 1e3
    return modis
