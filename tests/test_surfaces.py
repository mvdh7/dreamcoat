import pandas as pd
import numpy as np
import dreamcoat as dc
from cartopy import crs as ccrs

# Prepare adjacency edges for neutralocean surface calculation
edges = [
    [0, 24],
    [0, 27],
    [1, 8],
    [8, 19],
    [19, 24],
    [24, 27],
    [24, 35],
    [27, 35],
    [19, 40],
    [8, 46],
    [8, 50],
    [1, 59],
    [35, 40],
    [40, 46],
    [46, 50],
    [50, 59],
    [40, 70],
    [46, 63],
    [50, 74],
    [59, 112],
    [59, 117],
    [63, 70],
    [63, 74],
    [74, 112],
    [112, 117],
    [63, 93],
    [74, 89],
    [98, 117],
    [81, 86],
    [81, 112],
    [74, 81],
    [81, 89],
    [81, 98],
    [86, 89],
    [89, 93],
    [86, 98],
    [81, 117],
    [35, 70],
    [112, 126],
]
ctdz = pd.read_parquet("tests/data/ctdz.parquet")
graph = dc.CruiseGraph(ctdz, edges=edges)


def test_get_surfaces():
    station_ref = 1
    graph.get_surfaces(station_ref)
    assert station_ref in graph.surfaces_raw
    assert station_ref in graph.surfaces


def test_get_surfaces_all():
    graph.get_surfaces_all()
    for s in graph.stations.index:
        assert s in graph.surfaces


def test_plot_surface_map():
    extent = [-1, 7, 56.5, 64]
    crs = ccrs.EquidistantConic(central_longitude=np.mean([extent[0], extent[1]]))
    fig, ax = graph.plot_surface_map(1, 50, crs=crs, extent=extent)
    fig, ax = graph.plot_surface_map(8, 50, crs=crs, extent=extent)


# test_get_surfaces()
# test_get_surfaces_all()
# test_plot_surface_map()
