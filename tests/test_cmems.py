import numpy as np
import dreamcoat as dc
from matplotlib import pyplot as plt


# t8Ejkp-boxros-6qacgu-bUqmyr-i3jIls

lims = dict(
    date_min="2022-12-26",
    date_max="2022-12-30",
    longitude_min=-20,
    longitude_max=30,
    latitude_min=-50,
    latitude_max=-25,
)
# dc.cmems.download_surphys(
# filepath="tests/data/",
# **lims,
# )

surface = dc.cmems.open_surface(filepath="tests/data/", **lims)

speed_knots = 9
timestep_hours = 24
ship_distance = (
    dc.convert.knots_to_kph(speed_knots) * timestep_hours * np.array([1, 2, 3])
)

# surphys['dissic'] = surbio.dissic.interp_like(surphys.theta, method='nearest')

#%%
mooring_lon_lat = [10.02, -38.41]

for i in range(surface.time.size):
    fig, ax = dc.plot.surface_map(
        # surface.mean('time'),
        surface.isel(time=i),
        "ssh",
        land_visible=True,
        ship_lon_lat=mooring_lon_lat,
        ship_distance=ship_distance[0],
        # quiver_coarsen=20,
        # quiver_alpha=0.1,
        map_extent=[5, 15, -40, -35],
        # vmin=13,
        # vmax=24,
        # vmin=0,
        # vmax=0.3,
        # vmin=0,
        # vmax=80,
        # vmin=7.98,
        # vmax=8.14,
        # save_figure=True,
        # save_path="tests/figures/{}_".format(i),
    )
    plt.show()
    plt.close()

#%%
dc.plot.surface_map_daily(
    surface,
    "theta",
    save_figure=True,
    save_path="tests/figures/",
    # map_extent=[5, 15, -40, -35],
    # ship_lon_lat=mooring_lon_lat,
    ship_distance=ship_distance[0],
)

#%%
fvar = "current"
f_surface = surface.sel(
    longitude=mooring_lon_lat[0], latitude=mooring_lon_lat[1], method="nearest"
)

fig, ax = dc.plot.surface_timeseries(f_surface, fvar)

# f_surface.sel(
#     longitude=mooring_lon_lat[0], latitude=mooring_lon_lat[1], method="nearest"
# ).isel(depth=0)[fvar].plot(ax=ax)

#%%
fig, axs = dc.plot.surface_timeseries_grid(f_surface)

#%%
plt.plot(f_surface.current_east, f_surface.current_north)
plt.axhline(0)
plt.axvline(0)

#%% Polar
dc.plot.surface_currents(f_surface)
