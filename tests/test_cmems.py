import numpy as np
import dreamcoat as dc
from matplotlib import pyplot as plt


# boxros-6qacgu-bUqmyr

lims = dict(
    date_min="2021-01-20",
    date_max="2021-02-10",
    longitude_min=-70,
    longitude_max=30,
    latitude_min=-40,
    latitude_max=-20,
)
# dc.cmems.download_surphys(
# filepath="tests/data/",
# **lims,
# )
surphys = dc.cmems.open_surphys(
    filepath="tests/data/",
    **lims,
)


speed_knots = 9
timestep_hours = 24
ship_distance = (
    dc.convert.knots_to_kph(speed_knots) * timestep_hours * np.array([1, 2, 3])
)

#%%
mooring_lon_lat = [10.02, -38.41]

for i in range(surphys.time.size):
    fig, ax = dc.plot.surphys_map(
        # surphys.mean('time'),
        surphys.isel(time=i),
        "mld",
        land_visible=False,
        ship_lon_lat=mooring_lon_lat,
        ship_distance=10,
        # quiver_coarsen=20,
        # quiver_alpha=0.1,
        map_extent=[5, 15, -40, -35],
        # vmin=13,
        # vmax=24,
        # vmin=0,
        # vmax=0.3,
        vmin=0,
        vmax=80,
    )
    plt.show()
    plt.close()


#%%

fvar = "ssh"
f_surphys = surphys.sel(
    longitude=mooring_lon_lat[0], latitude=mooring_lon_lat[1], method="nearest"
).isel(depth=0)

fig, ax = dc.plot.surphys_timeseries(f_surphys, fvar)

# f_surphys.sel(
#     longitude=mooring_lon_lat[0], latitude=mooring_lon_lat[1], method="nearest"
# ).isel(depth=0)[fvar].plot(ax=ax)

#%%
fig, axs = dc.plot.surphys_timeseries_grid(f_surphys)

#%%
plt.plot(f_surphys.current_east, f_surphys.current_north)
plt.axhline(0)
plt.axvline(0)

#%% Polar


dc.plot.surphys_currents(f_surphys)
