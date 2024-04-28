"""
dreamcoat.underway
==================
Read and process underway data files (specifically, from RV Pelagia).

Functions
---------
read_sql
    Import the data from an sql underway file and parse it into a usable format.
"""

import pandas as pd
import numpy as np
from matplotlib import dates as mdates


def read_sql(filename):
    """Import the data from an sql underway file and parse it into a usable format.

    Parameters
    ----------
    filename : str
        The sql file name (and path).

    Returns
    -------
    pd.DataFrame
        The underway dataset.
    """
    # Import and parse underway SQL file
    with open(filename, "r") as f:
        underway_raw = f.readlines()
    headers = []
    is_header = False
    underway = []
    for line in underway_raw:
        if is_header:
            try:
                headers.append(line.split("`")[1])
            except IndexError:
                is_header = False
        if line.startswith("CREATE TABLE `aggregated` ("):
            is_header = True
        data_line = "INSERT INTO `aggregated` VALUES ("
        if line.startswith(data_line):
            line_data = (
                line.replace("'", "").split(data_line)[1].split(");\n")[0].split("),(")
            )
            line_data = [d.split(",") for d in line_data]
            underway = [*underway, *line_data]
    underway_array = np.array(underway)
    underway = pd.DataFrame({h: underway_array[:, i] for i, h in enumerate(headers)})
    # Convert all float columns to floats
    underway = underway.replace(["INVALID", ""], "-999")
    underway = underway.replace("°", "")
    for c in underway.columns:
        try:
            underway[c] = underway[c].astype(float)
            underway.loc[underway[c] == -999, c] = np.nan
        except ValueError:
            pass
    # Pythonise column names
    underway.columns = [
        c.lower().replace(" ", "_").replace("(", "_").replace(")", "")
        for c in underway.columns
    ]
    # Reformat where required
    underway["datetime"] = pd.to_datetime(underway.date)
    underway["datenum"] = mdates.date2num(underway.datetime)
    underway["latitude"] = underway.latitude.where(
        underway.latitude_direction == "N", -underway.latitude
    )
    underway["longitude"] = underway.longitude.where(
        underway.longitude_direction == "E", -underway.longitude
    )
    # Drop unnecessary columns
    underway = underway.drop(
        columns=[
            "date",
            "latitude_direction",
            "longitude_direction",
            "year",
            "month",
            "day",
        ]
    )
    return underway
