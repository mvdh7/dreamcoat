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

Modules
-------
convert
    Convert between different units and formats.
ctd
    Dealing with CTD data files.
glider
    Working with Slocum gliders.
maps
    Working with data in (longitude, latitude) space.
meta
    Metadata for dreamcoat.
plot
    Shortcuts for data visualisations.
underway
    Working with underway data files.

Classes
-------
LatLon
    A tool to convert between decimal degrees and degrees decimal minutes.
Route
    A set of waypoints defining a route.

Functions
---------
hello
    Report the version number.
"""

from . import convert, ctd, glider, maps, meta, plot, stats, underway
from .maps import Route
from .maps.degrees_decimal_minutes import LatLon
from .meta import __version__, hello
