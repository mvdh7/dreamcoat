import os
from datetime import date
import numpy as np
import xarray as xr
import gsw
import calkulate as calk
import koolstof as ks
from . import convert, meta

rename_phys = {
    "mlotst": "mld",
    "so": "salinity",
    "thetao": "theta",
    "uo": "current_east",
    "vo": "current_north",
    "zos": "ssh",
}
rename_bio = {
    "talk": "alkalinity",
    "dissic": "dic",
    "ph": "pH",
    "no3": "nitrate",
    "po4": "phosphate",
    "si": "silicate",
    "fe": "iron",
    "chl": "chlorophyll",
    "o2": "oxygen",
    "spco2": "pCO2",
    "nppv": "npp",
    "phyc": "phytoplankton",
}


def _get_filename_ending(
    date_min=None,
    date_max=None,
    longitude_min=None,
    longitude_max=None,
    latitude_min=None,
    latitude_max=None,
    extension="zarr",
):
    """Generate final, dataset-independent part of filename used to save downloaded
    CMEMS files.

    Parameters
    ----------
    date_min : str
        First date of data to download in '%Y-%m-%d' format.
    date_max : str
        Last date of data to download in '%Y-%m-%d' format.
    latitude_min : int
        Minimum latitude in decimal degrees N.
    latitude_max : int
        Maximum latitude in decimal degrees N.
    longitude_min : int
        Minimum longitude in decimal degrees E.
    longitude_max : int
        Maximum longitude in decimal degrees E.
    extension : str, optional
        The file extension, by default "zarr".

    Returns
    -------
    str
        The filename for the specified date and location ranges.
    """
    return "{}_{}_{}_{}_{}_{}.{}".format(
        date_min,
        date_max,
        longitude_min,
        longitude_max,
        latitude_min,
        latitude_max,
        extension,
    )


def _get_surphys_filename(
    date_min,
    date_max,
    longitude_min,
    longitude_max,
    latitude_min,
    latitude_max,
    extension="zarr",
):
    """Generate filename used to save downloaded CMEMS files from the
    GLOBAL_ANALYSIS_FORECAST_PHY_001_024 dataset.

    Parameters
    ----------
    date_min : str
        First date of data to download in '%Y-%m-%d' format.  Uses today if None.
    date_max : str
        Last date of data to download in '%Y-%m-%d' format.  Uses today if None.
    latitude_min : int
        Minimum latitude in decimal degrees N.
    latitude_max : int
        Maximum latitude in decimal degrees N.
    longitude_min : int
        Minimum longitude in decimal degrees E.
    longitude_max : int
        Maximum longitude in decimal degrees E.
    extension : str, optional
        The file extension, by default "zarr".

    Returns
    -------
    str
        The filename for the specified date and location ranges.
    """
    if date_min is None:
        date_min = date.today().strftime("%Y-%m-%d")
    if date_max is None:
        date_max = date.today().strftime("%Y-%m-%d")
    return "global-analysis-forecast-phy-001-024_" + _get_filename_ending(
        date_min,
        date_max,
        longitude_min,
        longitude_max,
        latitude_min,
        latitude_max,
        extension,
    )


def _get_deepphys_filename(
    date_min,
    date_max,
    longitude_min,
    longitude_max,
    latitude_min,
    latitude_max,
    extension="zarr",
):
    """Generate filename used to save 3D downloaded CMEMS files from the
    GLOBAL_ANALYSIS_FORECAST_PHY_001_024 dataset.

    Parameters
    ----------
    date_min : str
        First date of data to download in '%Y-%m-%d' format.  Uses today if None.
    date_max : str
        Last date of data to download in '%Y-%m-%d' format.  Uses today if None.
    latitude_min : int
        Minimum latitude in decimal degrees N.
    latitude_max : int
        Maximum latitude in decimal degrees N.
    longitude_min : int
        Minimum longitude in decimal degrees E.
    longitude_max : int
        Maximum longitude in decimal degrees E.
    extension : str, optional
        The file extension, by default "zarr".

    Returns
    -------
    str
        The filename for the specified date and location ranges.
    """
    if date_min is None:
        date_min = date.today().strftime("%Y-%m-%d")
    if date_max is None:
        date_max = date.today().strftime("%Y-%m-%d")
    return "global-analysis-forecast-phy-001-024_3D_" + _get_filename_ending(
        date_min,
        date_max,
        longitude_min,
        longitude_max,
        latitude_min,
        latitude_max,
        extension,
    )


def _get_surbio_filename(
    date_min,
    date_max,
    longitude_min,
    longitude_max,
    latitude_min,
    latitude_max,
    extension="zarr",
):
    """Generate filename used to save downloaded CMEMS files from the
    GLOBAL_ANALYSIS_FORECAST_BIO_001_028 dataset.

    Parameters
    ----------
    date_min : str
        First date of data to download in '%Y-%m-%d' format.  Uses today if None.
    date_max : str
        Last date of data to download in '%Y-%m-%d' format.  Uses today if None.
    latitude_min : int
        Minimum latitude in decimal degrees N.
    latitude_max : int
        Maximum latitude in decimal degrees N.
    longitude_min : int
        Minimum longitude in decimal degrees E.
    longitude_max : int
        Maximum longitude in decimal degrees E.
    extension : str, optional
        The file extension, by default "zarr".

    Returns
    -------
    str
        The filename for the specified date and location ranges.
    """
    if date_min is None:
        date_min = date.today().strftime("%Y-%m-%d")
    if date_max is None:
        date_max = date.today().strftime("%Y-%m-%d")
    return "global-analysis-forecast-bio-001-028_" + _get_filename_ending(
        date_min,
        date_max,
        longitude_min,
        longitude_max,
        latitude_min,
        latitude_max,
        extension,
    )


def _get_deepbio_filename(
    date_min,
    date_max,
    longitude_min,
    longitude_max,
    latitude_min,
    latitude_max,
    extension="zarr",
):
    """Generate filename used to save 3D downloaded CMEMS files from the
    GLOBAL_ANALYSIS_FORECAST_BIO_001_028 dataset.

    Parameters
    ----------
    date_min : str
        First date of data to download in '%Y-%m-%d' format.  Uses today if None.
    date_max : str
        Last date of data to download in '%Y-%m-%d' format.  Uses today if None.
    latitude_min : int
        Minimum latitude in decimal degrees N.
    latitude_max : int
        Maximum latitude in decimal degrees N.
    longitude_min : int
        Minimum longitude in decimal degrees E.
    longitude_max : int
        Maximum longitude in decimal degrees E.
    extension : str, optional
        The file extension, by default "zarr".

    Returns
    -------
    str
        The filename for the specified date and location ranges.
    """
    if date_min is None:
        date_min = date.today().strftime("%Y-%m-%d")
    if date_max is None:
        date_max = date.today().strftime("%Y-%m-%d")
    return "global-analysis-forecast-bio-001-028_3D_" + _get_filename_ending(
        date_min,
        date_max,
        longitude_min,
        longitude_max,
        latitude_min,
        latitude_max,
        extension,
    )


def download_surphys(
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
    convert_nc=True,
    delete_nc=True,
    get_current_vertical=False,
):
    """Download a CMEMS data file from the GLOBAL_ANALYSIS_FORECAST_PHY_001_024
    dataset including surface fields of salinity (so), potential temperature (thetao),
    mixed layer depth (mlotst), sea surface height (zos), and eastwards (uo) and
    northwards (vo) current velocities.

    Parameters
    ----------
    filename : str, optional
        File name to save to, by default None, in which case it is generated by
        _get_surphys_filename().
    filepath : str, optional
        File path for the saved file, NOT ending with a separator, by default "".
    date_min : str, optional
        First date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    date_max : str, optional
        Last date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    latitude_min : int, optional
        Minimum latitude in decimal degrees N, by default -90.
    latitude_max : int, optional
        Maximum latitude in decimal degrees N, by default 90.
    longitude_min : int, optional
        Minimum longitude in decimal degrees E, by default -180.
    longitude_max : int, optional
        Maximum longitude in decimal degrees E, by default 180.
    username : str, optional
        Your CMEMS username, by default None, in which case you will be prompted for it.
    password : str, optional
        Your CMEMS password, by default None, in which case you will be prompted for it.
    convert_nc : bool, optional
        Whether to convert the nc files into zarr/zip, by default True.
    delete_nc : bool, optional
        Whether to delete the nc files after converting, by default True.
    get_current_vertical : bool, optional
        Whether to also download the vertical current velocity, by default False.
    """
    # Deal with None inputs
    if username is None:
        username = meta.get_dat_data("cmems_username")
    if password is None:
        password = meta.get_dat_data("cmems_password")
    if date_min is None:
        date_min = date.today().strftime("%Y-%m-%d")
    if date_max is None:
        date_max = date.today().strftime("%Y-%m-%d")
    if filename is None:
        filename = _get_surphys_filename(
            date_min,
            date_max,
            longitude_min,
            longitude_max,
            latitude_min,
            latitude_max,
            extension="nc",
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
    # Now get vertical current velocity too, if requested
    if get_current_vertical:
        os.system(
            "motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu "
            + "--service-id GLOBAL_ANALYSISFORECAST_PHY_001_024-TDS "
            + "--product-id cmems_mod_glo_phy-wcur_anfc_0.083deg_P1D-m "
            + "--longitude-min {} ".format(longitude_min)
            + "--longitude-max {} ".format(longitude_max)
            + "--latitude-min {} ".format(latitude_min)
            + "--latitude-max {} ".format(latitude_max)
            + '--date-min "{} 12:00:00" '.format(date_min)
            + '--date-max "{} 12:00:00" '.format(date_max)
            + "--depth-min 0.49402499198913574 "
            + "--depth-max 0.49402499198913574 "
            + "--variable wo "
            + "--out-dir {} ".format(filepath)
            + "--out-name {} ".format(filename.replace(".nc", "_wo.nc"))
            + "--user {} --pwd {}".format(username, password)
        )
    # Convert from nc to zarr and a zip archive
    if convert_nc:
        convert.nc_to_zarr_zip(filepath + os.sep + filename, delete_nc=delete_nc)
        if get_current_vertical:
            convert.nc_to_zarr_zip(
                filepath + os.sep + filename.replace(".nc", "_wo.nc"),
                delete_nc=delete_nc,
            )


def download_deepphys(
    filename=None,
    filepath="",
    date_min=None,
    date_max=None,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    depth_min=0.1,
    depth_max=10000,
    username=None,
    password=None,
    convert_nc=True,
    delete_nc=True,
):
    """Download a CMEMS data file from the GLOBAL_ANALYSIS_FORECAST_PHY_001_024
    dataset including 3D fields of salinity (so), potential temperature (thetao),
    mixed layer depth (mlotst), sea surface height (zos), and eastwards (uo) and
    northwards (vo) current velocities.

    Parameters
    ----------
    filename : str, optional
        File name to save to, by default None, in which case it is generated by
        _get_deepphys_filename().
    filepath : str, optional
        File path for the saved file, NOT ending with a separator, by default "".
    date_min : str, optional
        First date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    date_max : str, optional
        Last date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    latitude_min : int, optional
        Minimum latitude in decimal degrees N, by default -90.
    latitude_max : int, optional
        Maximum latitude in decimal degrees N, by default 90.
    longitude_min : int, optional
        Minimum longitude in decimal degrees E, by default -180.
    longitude_max : int, optional
        Maximum longitude in decimal degrees E, by default 180.
    username : str, optional
        Your CMEMS username, by default None, in which case you will be prompted for it.
    password : str, optional
        Your CMEMS password, by default None, in which case you will be prompted for it.
    convert_nc : bool, optional
        Whether to convert the nc files into zarr/zip, by default True.
    delete_nc : bool, optional
        Whether to delete the nc files after converting, by default True.
    """
    # Deal with None inputs
    if username is None:
        username = meta.get_dat_data("cmems_username")
    if password is None:
        password = meta.get_dat_data("cmems_password")
    if date_min is None:
        date_min = date.today().strftime("%Y-%m-%d")
    if date_max is None:
        date_max = date.today().strftime("%Y-%m-%d")
    if filename is None:
        filename = _get_deepphys_filename(
            date_min,
            date_max,
            longitude_min,
            longitude_max,
            latitude_min,
            latitude_max,
            extension="nc",
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
        + "--depth-min {} --depth-max {} ".format(depth_min, depth_max)
        + "--variable mlotst --variable so --variable thetao "
        + "--variable uo --variable vo --variable zos "
        + "--out-dir {} ".format(filepath)
        + "--out-name {} ".format(filename)
        + "--user {} --pwd {}".format(username, password)
    )
    # Convert from nc to zarr and a zip archive
    if convert_nc:
        convert.nc_to_zarr_zip(filepath + os.sep + filename, delete_nc=delete_nc)


def download_surbio(
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
    convert_nc=True,
    delete_nc=True,
):
    """Download a missing CMEMS data file from the GLOBAL_ANALYSIS_FORECAST_BIO_001_028
    # dataset including surface fields of salinity (so), potential temperature (thetao),
    # mixed layer depth (mlotst), sea surface height (zos), and eastwards (uo) and
    # northwards (vo) current velocities.

    Parameters
    ----------
    filename : str, optional
        File name to save to, by default None, in which case it is generated by
        _get_surbio_filename().
    filepath : str, optional
        File path for the saved file, NOT ending with a separator, by default "".
    date_min : str, optional
        First date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    date_max : str, optional
        Last date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    latitude_min : int, optional
        Minimum latitude in decimal degrees N, by default -90.
    latitude_max : int, optional
        Maximum latitude in decimal degrees N, by default 90.
    longitude_min : int, optional
        Minimum longitude in decimal degrees E, by default -180.
    longitude_max : int, optional
        Maximum longitude in decimal degrees E, by default 180.
    username : str, optional
        Your CMEMS username, by default None, in which case you will be prompted for it.
    password : str, optional
        Your CMEMS password, by default None, in which case you will be prompted for it.
    convert_nc : bool, optional
        Whether to convert the nc file into a zarr/zip, by default True.
    """
    # Deal with None inputs
    if username is None:
        username = meta.get_dat_data("cmems_username")
    if password is None:
        password = meta.get_dat_data("cmems_password")
    if date_min is None:
        date_min = date.today().strftime("%Y-%m-%d")
    if date_max is None:
        date_max = date.today().strftime("%Y-%m-%d")
    if filename is None:
        filename = _get_surbio_filename(
            date_min,
            date_max,
            longitude_min,
            longitude_max,
            latitude_min,
            latitude_max,
            extension="nc",
        )
    # Download the file
    os.system(
        "motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu "
        + "--service-id GLOBAL_ANALYSIS_FORECAST_BIO_001_028-TDS "
        + "--product-id global-analysis-forecast-bio-001-028-daily "
        + "--longitude-min {} ".format(longitude_min)
        + "--longitude-max {} ".format(longitude_max)
        + "--latitude-min {} ".format(latitude_min)
        + "--latitude-max {} ".format(latitude_max)
        + '--date-min "{} 12:00:00" '.format(date_min)
        + '--date-max "{} 12:00:00" '.format(date_max)
        + "--depth-min 0.494 --depth-max 0.4941 "
        + "--variable chl --variable dissic --variable fe "
        + "--variable no3 --variable nppv --variable o2 "
        + "--variable ph --variable phyc --variable po4 "
        + "--variable si --variable spco2 --variable talk "
        + "--out-dir {} ".format(filepath)
        + "--out-name {} ".format(filename)
        + "--user {} --pwd {}".format(username, password)
    )
    # Convert from nc to zarr and a zip archive
    if convert_nc:
        convert.nc_to_zarr_zip(filepath + os.sep + filename, delete_nc=delete_nc)


def download_deepbio(
    filename=None,
    filepath="",
    date_min=None,
    date_max=None,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    depth_min=0.1,
    depth_max=10000,
    username=None,
    password=None,
    convert_nc=True,
    delete_nc=True,
):
    """Download a missing CMEMS data file from the GLOBAL_ANALYSIS_FORECAST_BIO_001_028
    # dataset including surface fields of salinity (so), potential temperature (thetao),
    # mixed layer depth (mlotst), sea surface height (zos), and eastwards (uo) and
    # northwards (vo) current velocities.

    Parameters
    ----------
    filename : str, optional
        File name to save to, by default None, in which case it is generated by
        _get_deepbio_filename().
    filepath : str, optional
        File path for the saved file, NOT ending with a separator, by default "".
    date_min : str, optional
        First date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    date_max : str, optional
        Last date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    latitude_min : int, optional
        Minimum latitude in decimal degrees N, by default -90.
    latitude_max : int, optional
        Maximum latitude in decimal degrees N, by default 90.
    longitude_min : int, optional
        Minimum longitude in decimal degrees E, by default -180.
    longitude_max : int, optional
        Maximum longitude in decimal degrees E, by default 180.
    username : str, optional
        Your CMEMS username, by default None, in which case you will be prompted for it.
    password : str, optional
        Your CMEMS password, by default None, in which case you will be prompted for it.
    convert_nc : bool, optional
        Whether to convert the nc file into a zarr/zip, by default True.
    """
    # Deal with None inputs
    if username is None:
        username = meta.get_dat_data("cmems_username")
    if password is None:
        password = meta.get_dat_data("cmems_password")
    if date_min is None:
        date_min = date.today().strftime("%Y-%m-%d")
    if date_max is None:
        date_max = date.today().strftime("%Y-%m-%d")
    if filename is None:
        filename = _get_deepbio_filename(
            date_min,
            date_max,
            longitude_min,
            longitude_max,
            latitude_min,
            latitude_max,
            extension="nc",
        )
    # Download the file
    os.system(
        "motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu "
        + "--service-id GLOBAL_ANALYSIS_FORECAST_BIO_001_028-TDS "
        + "--product-id global-analysis-forecast-bio-001-028-daily "
        + "--longitude-min {} ".format(longitude_min)
        + "--longitude-max {} ".format(longitude_max)
        + "--latitude-min {} ".format(latitude_min)
        + "--latitude-max {} ".format(latitude_max)
        + '--date-min "{} 12:00:00" '.format(date_min)
        + '--date-max "{} 12:00:00" '.format(date_max)
        + "--depth-min {} --depth-max {} ".format(depth_min, depth_max)
        + "--variable chl --variable dissic --variable fe "
        + "--variable no3 --variable nppv --variable o2 "
        + "--variable ph --variable phyc --variable po4 "
        + "--variable si --variable spco2 --variable talk "
        + "--out-dir {} ".format(filepath)
        + "--out-name {} ".format(filename)
        + "--user {} --pwd {}".format(username, password)
    )
    # Convert from nc to zarr and a zip archive
    if convert_nc:
        convert.nc_to_zarr_zip(filepath + os.sep + filename, delete_nc=delete_nc)


def open_surphys(
    filepath="",
    date_min=None,
    date_max=None,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    with_current_vertical=True,
):
    """Open a CMEMS data file from the GLOBAL_ANALYSIS_FORECAST_PHY_001_024 dataset.
    Output dataset includes surface fields of practical salinity (salinity), potential
    temperature in °C (theta), mixed layer depth in m (mld), sea surface height in m
    (ssh), eastwards (current_east) and northwards (current_north) current velocities in
    m/s, total horizontal current speed in m/s (current_speed), and vertical current
    speed in m/s (current_vertical).

    Parameters
    ----------
    filepath : str, optional
        File path for the saved file, by default "".
    date_min : str, optional
        First date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    date_max : str, optional
        Last date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    latitude_min : int, optional
        Minimum latitude in decimal degrees N, by default -90.
    latitude_max : int, optional
        Maximum latitude in decimal degrees N, by default 90.
    longitude_min : int, optional
        Minimum longitude in decimal degrees E, by default -180.
    longitude_max : int, optional
        Maximum longitude in decimal degrees E, by default 180.
    with_current_vertical : bool, optional
        Whether to import vertical current speed, by default True.

    Returns
    -------
    cmems
    """
    filename = _get_surphys_filename(
        date_min, date_max, longitude_min, longitude_max, latitude_min, latitude_max
    )
    cmems = xr.open_dataset(filepath + filename, engine="zarr").isel(depth=0)
    cmems["current_speed"] = np.sqrt(cmems.uo**2 + cmems.vo**2)
    cmems = cmems.rename(rename_cmems_phys)
    # Add the vertical current speed if requested, which has to be downloaded separately
    if with_current_vertical:
        cmems["current_vertical"] = (
            xr.open_dataset(
                filepath + filename.replace(".zarr", "_wo.zarr"), engine="zarr"
            )
            .wo.isel(depth=0)
            .interp_like(cmems, method="nearest")
        )
    return cmems


def open_surbio(
    filepath="",
    date_min=None,
    date_max=None,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
):
    """Open a CMEMS data file from the GLOBAL_ANALYSIS_FORECAST_BIO_001_028 dataset.

    Parameters
    ----------
    filepath : str, optional
        File path for the saved file, by default "".
    date_min : str, optional
        First date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    date_max : str, optional
        Last date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    latitude_min : int, optional
        Minimum latitude in decimal degrees N, by default -90.
    latitude_max : int, optional
        Maximum latitude in decimal degrees N, by default 90.
    longitude_min : int, optional
        Minimum longitude in decimal degrees E, by default -180.
    longitude_max : int, optional
        Maximum longitude in decimal degrees E, by default 180.

    Returns
    -------
    cmems
    """
    filename = _get_surbio_filename(
        date_min, date_max, longitude_min, longitude_max, latitude_min, latitude_max
    )
    return xr.open_dataset(filepath + filename, engine="zarr").isel(depth=0)


def open_surface(
    filepath="",
    date_min=None,
    date_max=None,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    with_current_vertical=False,
):
    """Open CMEMS data files from the GLOBAL_ANALYSIS_FORECAST_PHY_001_024 and
    GLOBAL_ANALYSIS_FORECAST_BIO_001_028 datasets.  Combine both datasets together,
    nearest-neighbour interpolating the lower-resolution biogeochemical dataset to
    match the physical dataset's spatial grid.  Adjust the units of the biogeochemical
    variables to more usual ones.

    Parameters
    ----------
    filepath : str, optional
        File path for the saved files, by default "".
    date_min : str, optional
        First date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    date_max : str, optional
        Last date of data to download in '%Y-%m-%d' format, by default None, in which
        case today is used.
    latitude_min : int, optional
        Minimum latitude in decimal degrees N, by default -90.
    latitude_max : int, optional
        Maximum latitude in decimal degrees N, by default 90.
    longitude_min : int, optional
        Minimum longitude in decimal degrees E, by default -180.
    longitude_max : int, optional
        Maximum longitude in decimal degrees E, by default 180.
    with_current_vertical : bool, optional
        Whether to import vertical current speed, by default False.

    Returns
    -------
    surface
        An xarray Dataset containing the requested data.
    """
    # Open the physical and biological datasets
    surface = open_surphys(
        filepath=filepath,
        date_min=date_min,
        date_max=date_max,
        latitude_min=latitude_min,
        latitude_max=latitude_max,
        longitude_min=longitude_min,
        longitude_max=longitude_max,
        with_current_vertical=with_current_vertical,
    )
    surface_bio = open_surbio(
        filepath=filepath,
        date_min=date_min,
        date_max=date_max,
        latitude_min=latitude_min,
        latitude_max=latitude_max,
        longitude_min=longitude_min,
        longitude_max=longitude_max,
    )
    # Copy all surface_bio fields into surface, interpolating with nearest neighbour
    # to match the higher resolution of the physical dataset
    for k, v in surface_bio.items():
        surface[k] = v.interp_like(surface.theta, method="nearest")
    # Convert concentrations (per m3) to substance contents (per kg)
    surface["density"] = calk.density.seawater_1atm_MP81(
        temperature=surface.theta, salinity=surface.salinity
    )
    surface["density"].attrs.update(
        {
            "units": "kg dm-3",
            "unit_long": "kilograms per cubic decimetre",
        }
    )
    surface["density_anomaly"] = (surface.density - 1) * 1e3
    dvars = [
        "o2",
        "no3",
        "po4",
        "si",
        "talk",
        "dissic",
        "fe",
    ]
    for v in dvars:
        surface[v] /= surface.density
    # Adjust other units to PyCO2SYS standards
    surface["talk"] *= 1e3  # now µmol/kg
    surface["dissic"] *= 1e3  # now µmol/kg
    surface["spco2"] /= 101325e-6  # now µatm
    surface["fe"] *= 1e3  # now nmol/kg
    # Update Dataset attributes to reflect adjusted units
    umolkg = ["talk", "dissic", "o2", "no3", "po4", "si"]
    for v in umolkg:
        surface[v].attrs.update(
            {
                "units": "µmol kg-1",
                "unit_long": "micromoles per kilogram seawater",
            }
        )
    surface["fe"].attrs.update(
        {
            "units": "nmol kg-1",
            "unit_long": "nanomoles per kilogram seawater",
        }
    )
    surface["spco2"].attrs.update(
        {
            "units": "µatm",
            "unit_long": "microatmospheres",
        }
    )
    # Calculate additional variables of interest
    surface["aou"] = ks.aou_GG92(
        oxygen=surface.o2, temperature=surface.theta, salinity=surface.salinity
    )[0]
    surface["aou"].attrs.update(
        {
            "long_name": "Apparent oxygen utilisation",
            "units": "µmol kg-1",
            "unit_long": "micromoles per kilogram seawater",
        }
    )
    return surface


def open_phys_and_bio(
    filename_phys="cmems_phys.zarr", filename_bio="cmems_bio.zarr", filepath="."
):
    """Read a pair of CMEMS physical and biogeochemical files that correspond to the
    same time and space domains together, merge into a single dataset and compute a
    bunch of extra variables of interest.

    Parameters
    ----------
    filename_phys : str, optional
        Filename for the CMEMS physical dataset, by default "cmems_phys.zarr".
    filename_bio : str, optional
        Filename for the CMEMS biogeochemical dataset, by default "cmems_bio.zarr".
    filepath : str, optional
        Path to both datasets excluding a trailing file separator, by default ".".

    Returns
    -------
    xarray.Dataset
        The combined CMEMS dataset.
    """
    # Import CMEMS datasets
    cmems = xr.open_dataset(
        os.sep.join((filepath, filename_phys)), engine="zarr"
    ).rename(rename_phys)
    cmems_bio = xr.open_dataset(
        os.sep.join((filepath, filename_bio)),
        engine="zarr",
    ).rename(rename_bio)
    # Do GSW calculations
    cmems["pressure"] = gsw.p_from_z(-cmems.depth, cmems.latitude)
    cmems["pressure"].attrs.update(
        {
            "long_name": "Hydrostatic pressure",
            "standard_name": "pressure",
            "unit_long": "decibar",
            "units": "dbar",
        }
    )
    cmems["salinity_absolute"] = gsw.SA_from_SP(
        cmems.salinity, cmems.pressure, cmems.longitude, cmems.latitude
    )
    cmems["salinity_absolute"].attrs.update(
        {
            "long_name": "Absolute salinity",
            "standard_name": "salinity_absolute",
            "unit_long": "grams per kilogram",
            "units": "g/kg",
        }
    )
    cmems["temperature_conservative"] = gsw.CT_from_pt(
        cmems.salinity_absolute, cmems.theta
    )
    cmems["temperature_conservative"].attrs.update(
        {
            "long_name": "Conservative temperature",
            "standard_name": "temperature_conservative",
            "unit_long": "degrees Celsius",
            "units": "degrees_C",
        }
    )
    cmems["density"] = (
        gsw.rho(cmems.salinity_absolute, cmems.temperature_conservative, cmems.pressure)
        / 1000
    )
    cmems["density"].attrs.update(
        {
            "long_name": "Density in situ",
            "standard_name": "density",
            "unit_long": "kilograms per litre",
            "units": "kg/L",
        }
    )
    # Move bio variables across to main dataset and interpolate them onto the same grid
    for k, v in cmems_bio.items():
        cmems[k] = v.interp_like(cmems.theta, method="nearest")
    # Convert concentrations (per m3) to substance contents (per kg)
    dvars = [
        "oxygen",
        "nitrate",
        "phosphate",
        "silicate",
        "alkalinity",
        "dic",
        "iron",
    ]
    for v in dvars:
        cmems[v] /= cmems.density
    # Adjust other units to PyCO2SYS standards
    cmems["alkalinity"] *= 1e3  # now µmol/kg
    cmems["dic"] *= 1e3  # now µmol/kg
    cmems["pCO2"] /= 101325e-6  # now µatm
    cmems["iron"] *= 1e3  # now nmol/kg
    # Update Dataset attributes to reflect adjusted units
    umolkg = ["alkalinity", "dic", "oxygen", "nitrate", "phosphate", "silicate"]
    for v in umolkg:
        cmems[v].attrs.update(
            {
                "units": "µmol kg-1",
                "unit_long": "micromoles per kilogram seawater",
            }
        )
    cmems["iron"].attrs.update(
        {
            "units": "nmol kg-1",
            "unit_long": "nanomoles per kilogram seawater",
        }
    )
    cmems["pCO2"].attrs.update(
        {
            "units": "µatm",
            "unit_long": "microatmospheres",
        }
    )
    # Calculate additional variables of interest
    cmems["aou"] = ks.aou_GG92(
        oxygen=cmems.oxygen, temperature=cmems.theta, salinity=cmems.salinity
    )[0]
    cmems["aou"].attrs.update(
        {
            "long_name": "Apparent oxygen utilisation from potential temperature",
            "units": "µmol kg-1",
            "unit_long": "micromoles per kilogram seawater",
        }
    )
    return cmems
