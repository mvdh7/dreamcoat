"""Dreamcoat
=========

> The joyful oceanographic seagoing expedition planning helper

Dreamcoat is a set of tools to help work with live oceanographic data of the sort that
might be needed during a research expedition.  At sea, internet bandwidth is usually
restricted, so dreamcoat is designed to be remotely usable: it can be left to run on a
server on land, which does the heavy lifting in terms of downloading files, producing
smaller figures or data subsets that can be automatically sent by email or uploaded to
an online repo.

Dreamcoat is being developed by Dr Matthew Humphreys (https://seaco2.group) at the NIOZ
Royal Netherlands Institute for Sea Research, Texel.  You're very welcome to give it a
go and suggestions and contributions are welcome, but please be aware that it is
primarily designed as a personal tool so I probably won't be able to help out if you get
stuck; nothing is guaranteed not to break between versions, and the documentation
(https://mvdh.xyz/dreamcoat) may be incomplete or out of date.

The docstring examples assume that `dreamcoat` has been imported as ``dc``:

  >>> import dreamcoat as dc

Available subpackages
---------------------
cmems
    Tools to import and open CMEMS data files
    (https://data.marine.copernicus.eu/products).
convert
    Functions to convert between different units and formats.
ctd
    Dealing with CTD data files.
maps
    Working with data in (longitude, latitude) space.
meta
    Metadata for `dreamcoat`.
modis
    Tools to import and open data from the MODIS satellites
    (https://oceancolor.gsfc.nasa.gov/).
plot
    Shortcuts for data visualisations.
send
    Send files by email.
style
    Get colour schemes and nicely formatted labels based on variable names.
"""

from . import convert, cmems, ctd, maps, meta, modis, plot, send, style
from .maps import (
    Route,
    get_route_distance,
    linspace_gc,
    linspace_gc_waypoints,
    build_vptree,
)
from .maps.degrees_decimal_minutes import DegreesDecimalMinutes

DDM = DegreesDecimalMinutes
from .meta import __version__, hello
from .style import styles
