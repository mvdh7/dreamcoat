import pandas as pd
import numpy as np
from matplotlib import dates as mdates

ctd_renamer = {
    "scan": "scan",
    "prDM": "pressure",
    "t090C": "temperature",
    "CStarAt0": "beam_attenuation",
    "CStarTr0": "beam_transmission",
    "flC": "fluorescence_chelsea",
    "v2": "voltage_2",
    "ph": "pH_voltage",
    "sbeox0Mg/L": "oxygen_mgl",
    "sbeox0Mm/L": "oxygen_mml",
    "altM": "altimeter",
    "flECO-AFL": "fluorescence_wetlabs",
    "turbWETntu0": "turbidity",
    "timeS": "elapsed_seconds",
    "c0mS/cm": "conductivity",
    "sal00": "salinity",
    "potemp090C": "theta",
    "sigma-é00": "density",
    "depSM": "depth",
    "svCM": "sound_velocity",
    "bottle": "bottle",
    "c0S/m": "conductivity",
    "spar/sat/lin": "spar",
    "par/sat/log": "par",
    "sigma-t00": "sigma_theta",
}
ctd_Renamer = {k[0].upper() + k[1:]: v for k, v in ctd_renamer.items()}


def read_cnv_1Hz(filename):
    """Read the contents of a 1 Hz cnv file from a CTD cast.

    Parameters
    ----------

    filename : str
        The filename (and path) to read.

    Returns
    -------
    pd.DataFrame
        The imported dataset.
    """
    # Find number of header lines to skip and get column headers and other info
    with open(filename, "r", encoding="unicode_escape", errors="replace") as f:
        skiprows = 0
        names = []
        ctd_Hz_line = ""
        while not ctd_Hz_line.startswith("*END*"):
            ctd_Hz_line = f.readline()
            # Get column headers
            if ctd_Hz_line.startswith("# name "):
                names.append(ctd_Hz_line.split(" = ")[1].split(":")[0])
            # Get latitude and longitude
            if ctd_Hz_line.startswith("* NMEA Latitude = "):
                latitude_txt = ctd_Hz_line.split(" = ")[1].split(" ")
                latitude = float(latitude_txt[0]) + float(latitude_txt[1]) / 60
                if latitude_txt[2] == "S\n":
                    latitude *= -1
            if ctd_Hz_line.startswith("* NMEA Longitude = "):
                longitude_txt = ctd_Hz_line.split(" = ")[1].split(" ")
                longitude = float(longitude_txt[0]) + float(longitude_txt[1]) / 60
                if longitude_txt[2] == "W\n":
                    longitude *= -1
            # Get datetime information
            if ctd_Hz_line.startswith("# start_time = "):
                start_time = pd.to_datetime(ctd_Hz_line.split(" = ")[1].split("[")[0])
            # Increment row counter
            skiprows += 1
    # Deal with duplicate names by appending "_2" to their second appearance
    names_duplicate = []
    for i, name in enumerate(names):
        if names.count(name) == 2:
            names[names.index(name, i + 1)] += "_2"
            names_duplicate.append(names[names.index(name + "_2", i + 1)])
    # Import the dataset with pandas
    ctd_Hz = (
        pd.read_table(
            filename,
            encoding="unicode_escape",
            encoding_errors="replace",
            names=names,
            skiprows=skiprows,
            delim_whitespace=True,
        )
        .drop(columns=names_duplicate)
        .rename(columns=ctd_renamer)
    )
    # Add positional information to the dataframe
    ctd_Hz["datetime"] = start_time + np.arange(0, ctd_Hz.shape[0]).astype(
        "timedelta64[s]"
    )
    ctd_Hz["latitude"] = latitude
    ctd_Hz["longitude"] = longitude
    return ctd_Hz


def read_btl_raw(filename):
    # Find number of header lines to skip and get column headers and other info
    with open(filename, "r", encoding="ascii") as f:
        skiprows = -1
        names = []
        ctd_Hz_line = "*"
        while ctd_Hz_line.startswith("*") or ctd_Hz_line.startswith("#"):
            ctd_Hz_line = f.readline()
            skiprows += 1
    # Open the file
    btl_raw = (
        pd.read_fwf(filename, skiprows=skiprows, encoding="ascii")
        .rename(columns=ctd_Renamer)
        .iloc[1:]
    )
    # Parse datetime
    btl_raw["date"] = btl_raw.Date
    btl_raw["time"] = btl_raw.Date
    btl_raw.loc[btl_raw.bottle.isnull(), "date"] = np.nan
    btl_raw.loc[~btl_raw.bottle.isnull(), "time"] = np.nan
    btl_raw.date.fillna(method="pad", inplace=True)
    btl_raw.time.fillna(method="bfill", inplace=True)
    btl_raw["datetime"] = pd.to_datetime(btl_raw.date + "T" + btl_raw.time)
    btl_raw["datenum"] = mdates.date2num(btl_raw.datetime)
    btl_raw.drop(columns=["Date", "date", "time"], inplace=True)
    return btl_raw


def read_ctd_create_btl(filename_cnv_1Hz, filename_btl_raw):
    ctd = read_cnv_1Hz(filename_cnv_1Hz)
    btl_raw = read_btl_raw(filename_btl_raw)
    # Make new btl table
    btl = btl_raw[~btl_raw.bottle.isnull()][
        ["bottle", "datetime", "datenum"]
    ].set_index("bottle")
    # Add data from 30 seconds before firing for each sensor
    time_before_firing = pd.Timedelta("30s")
    for c in ctd.columns:
        if c not in btl:
            btl[c] = np.nan
            btl[c + "_std"] = np.nan
            btl[c + "_flag"] = np.nan
    ctd["bottle"] = 0
    for i, row in btl.iterrows():
        L = (ctd.datetime > row.datetime - time_before_firing) & (
            ctd.datetime <= row.datetime
        )
        ctd.loc[L, "bottle"] = i
        for c in ctd.columns:
            if c not in ["bottle", "datetime", "datenum"]:
                btl.loc[i, c] = ctd.loc[L, c].mean()
                btl.loc[i, c + "_std"] = ctd.loc[L, c].std()
    return ctd, btl
