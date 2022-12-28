import numpy as np
import dreamcoat as dc

# t8Ejkp-boxros-6qacgu-bUqmyr-i3jIls

lims = dict(
    date_min="2022-12-26",
    date_max="2022-12-30",
    longitude_min=-20,
    longitude_max=30,
    latitude_min=-50,
    latitude_max=-25,
)
surface = dc.cmems.open_surface(filepath="tests/data/", **lims)

speed_knots = 9
timestep_hours = 24
ship_distance = (
    dc.convert.knots_to_kph(speed_knots) * timestep_hours * np.array([1, 2, 3])
)
mooring_lon_lat = [10.02, -38.41]

#%% Daily surface maps
dc.plot.surface_map_daily(
    surface,
    "ph",
    save_figure=True,
    save_path="tests/figures/",
    # map_extent=[5, 15, -40, -35],
    # ship_lon_lat=mooring_lon_lat,
    ship_distance=ship_distance[0],
)

#%% Time-series plots
fvar = "dissic"
f_surface = surface.sel(
    longitude=mooring_lon_lat[0], latitude=mooring_lon_lat[1], method="nearest"
)
dc.plot.surface_timeseries_grids(f_surface)

#%% Currents at point as a polar plot
dc.plot.surface_currents(f_surface)
