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


class DegreesDecimalMinutes:
    def __init__(self, longitude=0, latitude=0):
        # Get values in decimal degrees
        if isinstance(longitude, str):
            self.longitude_dd = to_lon_dd(longitude)
        elif isinstance(longitude, float) or isinstance(longitude, int):
            self.longitude_dd = longitude
        if isinstance(latitude, str):
            self.latitude_dd = to_lat_dd(latitude)
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
        self.dd = self.longitude_dd, self.latitude_dd

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

    def to_dd(self):
        return self.longitude_dd, self.latitude_dd

    def __repr__(self):
        return "{:02.0f}°{:06.3f}'{}, {:03.0f}°{:06.3f}'{}".format(
            self.latitude_r,
            self.latitude_dm,
            self.latitude_dir,
            self.longitude_r,
            self.longitude_dm,
            self.longitude_dir,
        )

    def __add__(self, other):
        return DegreesDecimalMinutes(
            self.longitude_dd + other.longitude_dd, self.latitude_dd + other.latitude_dd
        )

    def __sub__(self, other):
        return DegreesDecimalMinutes(
            self.longitude_dd - other.longitude_dd, self.latitude_dd - other.latitude_dd
        )

    def __mul__(self, other):
        return DegreesDecimalMinutes(
            self.longitude_dd * other, self.latitude_dd * other
        )

    def __truediv__(self, other):
        return DegreesDecimalMinutes(
            self.longitude_dd / other, self.latitude_dd / other
        )
