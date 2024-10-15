import numpy as np


def _split_lon_ddm(lon_ddm):
    lon_split = lon_ddm.split("-")
    if len(lon_split) == 1:
        lon_split = [lon_ddm.split("°")[0], *lon_ddm.split("°")[1].split("'")]
    return lon_split


def _split_to_lon_dd(lon_split):
    return (float(lon_split[0]) + float(lon_split[1]) / 60) * dict(E=1, W=-1)[
        lon_split[2].upper()
    ]


def to_lon_dd(lon_ddm):
    return _split_to_lon_dd(_split_lon_ddm(lon_ddm))


def _split_lat_ddm(lat_ddm):
    lat_split = lat_ddm.split("-")
    if len(lat_split) == 1:
        lat_split = [lat_ddm.split("°")[0], *lat_ddm.split("°")[1].split("'")]
    return lat_split


def _split_to_lat_dd(lat_split):
    return (float(lat_split[0]) + float(lat_split[1]) / 60) * dict(N=1, S=-1)[
        lat_split[2].upper()
    ]


def to_lat_dd(lat_ddm):
    return _split_to_lat_dd(_split_lat_ddm(lat_ddm))


class LatLon:
    """A tool to convert between decimal degrees and degrees decimal minutes and to
    carry out basic calculations (e.g., add and subtract) using degrees decimal minutes.

    Parameters
    ----------
    latitude : float or str
        A latitude value.
          - If a float, should be in decimal degrees, with N positive.
          - If a str, should be in degrees decimal minutes, formatted either as (e.g.)
            "65°12.3'N" or "65-12.3-N".
    longitude : float or str
        A longitude value.
          - If a float, should be in decimal degrees, with E positive.
          - If a str, should be in degrees decimal minutes, formatted either as (e.g.)
            "65°12.3'E" or "65-12.3-E".

    Attributes
    ----------
    latitude_dd : float
        Latitude in decimal degrees, from -90 to 90 with N positive.
    latitude_r : int
        The floor of the latitude (degrees only).
    latitude_dm : float
        The decimal minutes of the latitude.
    latitude_dir : str
        The direction of the latitude ("N" or "S").
    longitude_dd : float
        Longitude in decimal degrees, from -180 to 180 with E positive.
    longitude_r : int
        The floor of the longitude (degrees only).
    longitude_dm : float
        The decimal minutes of the longitude.
    longitude_dir : str
        The direction of the longitude ("E" or "W").

    Methods
    -------
    to_dd : tuple of float
        Return (longitude, latitude) in decimal degrees.
    """

    def __init__(
        self,
        latitude=0,
        longitude=0,
    ):
        # Get values in decimal degrees
        if isinstance(longitude, str):
            self.longitude_dd = to_lon_dd(longitude.replace(" ", ""))
        elif isinstance(longitude, float) or isinstance(longitude, int):
            self.longitude_dd = longitude
        if isinstance(latitude, str):
            self.latitude_dd = to_lat_dd(latitude.replace(" ", ""))
        elif isinstance(latitude, float) or isinstance(latitude, int):
            self.latitude_dd = latitude
        # Adjust numbers to stay within valid ranges
        # (-90 to 90 for latitude, -180 to 180 for longitude)
        if self.latitude_dd > 90:
            self.latitude_dd = 180 - self.latitude_dd
            self.longitude_dd += 180
        self.keep_longitude_in_range()
        # Get rounded values and decimal minutes
        self.update_rounded_dm()

    def keep_longitude_in_range(self):
        if self.longitude_dd > 180:
            self.longitude_dd -= 360
        elif self.longitude_dd <= -180:
            self.longitude_dd += 360

    def update_rounded_dm(self):
        self.longitude_r = int(np.floor(np.abs(self.longitude_dd)))
        self.longitude_dm = 60 * (np.abs(self.longitude_dd) - self.longitude_r)
        self.longitude_dir = "EEW"[np.sign(self.longitude_dd).astype(int)]
        self.latitude_r = int(np.floor(np.abs(self.latitude_dd)))
        self.latitude_dm = 60 * (np.abs(self.latitude_dd) - self.latitude_r)
        self.latitude_dir = "NNS"[np.sign(self.latitude_dd).astype(int)]

    def to_dd(self, lon_first=True):
        """Return (longitude, latitude) in decimal degrees.

        Parameters
        ----------
        lon_first : bool, optional
            Whether to lead with longitude (True) or latitude (False), by default True.

        Returns
        -------
        tuple of float
            The (longitude, latitude) value in decimal degrees, with E and N positive
            (or (latitude, longitude) if not lon_first).
        """
        if lon_first:
            return self.longitude_dd, self.latitude_dd
        else:
            return self.latitude_dd, self.longitude_dd

    def to_ddm(self):
        """Return latitude and longitude in degrees decimal minutes.

        Returns
        -------
        str
            A nicely formatted string of the latitude and longitude values.
        """
        return "{:02.0f}°{:06.3f}'{}, {:03.0f}°{:06.3f}'{}".format(
            self.latitude_r,
            self.latitude_dm,
            self.latitude_dir,
            self.longitude_r,
            self.longitude_dm,
            self.longitude_dir,
        )

    def __repr__(self):
        return self.to_ddm() + " ({:.4f}, {:.4f})".format(
            self.longitude_dd, self.latitude_dd
        )

    def __add__(self, other):
        return LatLon(
            self.longitude_dd + other.longitude_dd, self.latitude_dd + other.latitude_dd
        )

    def __sub__(self, other):
        return LatLon(
            self.longitude_dd - other.longitude_dd, self.latitude_dd - other.latitude_dd
        )

    def __mul__(self, other):
        return LatLon(self.longitude_dd * other, self.latitude_dd * other)

    def __truediv__(self, other):
        return LatLon(self.longitude_dd / other, self.latitude_dd / other)


# def dd_to_ddm(dd):
#     """Convert decimal degrees into degrees decimal minutes.
#
#     Parameters
#     ----------
#     dd : float
#         (longitude, latitude) in decimal degrees.
#
#     Returns
#     -------
#     str
#         The position(s) in degrees decimal minutes.
#     """
#     single = dd.ndim == 1
#     if single:
#         dd = [dd]
#     ddm = []
#     for d in dd:
#         lon, lat = d
#         ddm.append(
#             "{:02.0f}°{:04.1f}'{}, {:03.0f}°{:04.1f}'{}".format(
#                 np.floor(np.abs(lat)),
#                 60 * (np.abs(lat) % 1),
#                 "SN"[int(lat > 0)],
#                 np.floor(np.abs(lon)),
#                 60 * (np.abs(lon) % 1),
#                 "WE"[int(lon > 0)],
#             )
#         )
#     if single:
#         ddm = ddm[0]
#     return ddm
