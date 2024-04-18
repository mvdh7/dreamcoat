"""
Dreamcoat
=========

> The joyful oceanographic seagoing expedition planning helper

Dreamcoat is a set of tools to help work with oceanographic data of the sort that might
be needed or generated during a research expedition.

Dreamcoat is being developed by Dr Matthew Humphreys (https://seaco2.group) at the NIOZ
Royal Netherlands Institute for Sea Research, Texel.  You're very welcome to give it a
go and suggestions and contributions are welcome, but please be aware that it is
primarily designed as a personal tool so there probably won't be much help if you get
stuck; nothing is guaranteed not to break between versions, and the documentation
(https://mvdh.xyz/dreamcoat) may be incomplete or out of date.

The docstring examples assume that `dreamcoat` has been imported as ``dc``:

  >>> import dreamcoat as dc

Available subpackages
---------------------
convert
    Convert between different units and formats.
ctd
    Dealing with CTD data files.
maps
    Working with data in (longitude, latitude) space.
meta
    Metadata for dreamcoat.
plot
    Shortcuts for data visualisations.
"""

from . import convert, ctd, maps, meta, plot
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
