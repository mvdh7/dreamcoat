import pandas as pd
import numpy as np


def read_goto(filename):
    """Import and parse route waypoints in a goto file.

    Parameters
    ----------
    filename : str
        The file name (and file path) to the goto file.

    Returns
    -------
    pd.DataFrame
        The route waypoints.
    """
    with open(filename, "r") as f:
        data = f.readlines()
    is_waypoints = False
    for i, l in enumerate(data):
        wpline = "b_arg: num_waypoints(nodim)"
        if l.lstrip().startswith("b_arg: num_waypoints(nodim)"):
            num_waypoints = int(l.lstrip().split(wpline)[1].split("#")[0])
        elif l.strip() == "<start:waypoints>":
            break
    waypoints = data[i + 1 : i + 1 + num_waypoints]
    waypoints = np.array([wp.split() for wp in waypoints])
    route = pd.DataFrame(
        {"longitude_text": waypoints[:, 0], "latitude_text": waypoints[:, 1]}
    )
    lon_spl = route.longitude_text.str.extract(r"([-]?)(\d+)(\d{2})\.(\d+)")
    route["longitude"] = np.array([1, -1])[(lon_spl[0] == "-").astype(int)] * (
        lon_spl[1].astype(float)
        + (lon_spl[2].astype(float) + lon_spl[3].astype(float) / 100) / 60
    )
    lat_spl = route.latitude_text.str.extract(r"([-]?)(\d+)(\d{2})\.(\d+)")
    route["latitude"] = np.array([1, -1])[(lat_spl[0] == "-").astype(int)] * (
        lat_spl[1].astype(float)
        + (lat_spl[2].astype(float) + lat_spl[3].astype(float) / 100) / 60
    )
    return route
