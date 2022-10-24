import os
from datetime import date
import xarray as xr


def get_filename(
    date_min, date_max, longitude_min, longitude_max, latitude_min, latitude_max
):
    if not date_min:
        date_min = date.today().strftime("%Y-%m-%d")
    if not date_max:
        date_max = date.today().strftime("%Y-%m-%d")
    return "global-analysis-forecast-phy-001-024_{}_{}_{}_{}_{}_{}.nc".format(
        date_min, date_max, longitude_min, longitude_max, latitude_min, latitude_max
    )


def download_cmems(
    filename=None,
    filepath="",
    date_min=None,
    date_max=None,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    username=None,
    password=None,
):
    # Deal with None inputs
    if not username:
        username = input("Please enter your CMEMS username: ")
    if not password:
        password = input("Please enter your CMEMS password: ")
    if not date_min:
        date_min = date.today().strftime("%Y-%m-%d")
    if not date_max:
        date_max = date.today().strftime("%Y-%m-%d")
    if not filename:
        filename = get_filename(
            date_min, date_max, longitude_min, longitude_max, latitude_min, latitude_max
        )
    # Download the file
    os.system(
        "motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu "
        + "--service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS "
        + "--product-id global-analysis-forecast-phy-001-024 "
        + "--longitude-min {} ".format(longitude_min)
        + "--longitude-max {} ".format(longitude_max)
        + "--latitude-min {} ".format(latitude_min)
        + "--latitude-max {} ".format(latitude_max)
        + '--date-min "{} 12:00:00" '.format(date_min)
        + '--date-max "{} 12:00:00" '.format(date_max)
        + "--depth-min 0.494 --depth-max 0.4941 "
        + "--variable mlotst --variable so --variable thetao "
        + "--variable uo --variable vo --variable zos "
        + "--out-dir {} ".format(filepath)
        + "--out-name {} ".format(filename)
        + "--user {} --pwd {}".format(username, password)
    )


def open_cmems(
    filepath="",
    date_min=None,
    date_max=None,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    username=None,
    password=None,
):
    filename = get_filename(
        date_min, date_max, longitude_min, longitude_max, latitude_min, latitude_max
    )
    try:
        cmems = xr.open_dataset(filepath + filename)
    except FileNotFoundError:
        download_cmems(
            filename=filename,
            filepath=filepath,
            date_min=date_min,
            date_max=date_max,
            latitude_min=latitude_min,
            latitude_max=latitude_max,
            longitude_min=longitude_min,
            longitude_max=longitude_max,
            username=username,
            password=password,
        )
        cmems = xr.open_dataset(filepath + filename)
    return cmems
