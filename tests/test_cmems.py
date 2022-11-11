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
        "theta",
        land_visible=False,
        ship_lon_lat=mooring_lon_lat,
        ship_distance=10,
        # quiver_coarsen=20,
        # quiver_alpha=0.1,
        map_extent=[5, 15, -40, -35],
        vmin=13,
        vmax=24,
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
plt.scatter(f_surphys.current_east, f_surphys.current_north)
plt.axhline(0)
plt.axvline(0)

#%% Polar
# from matplotlib import dates as mdates, pyplot as plt, patheffects as pe
# from matplotlib.collections import LineCollection
# from scipy.interpolate import pchip_interpolate
# from dreamcoat import convert

dc.plot.surphys_currents(f_surphys)


# def _get_surphys_currents_line(data, interpolate_pchip):
#     if interpolate_pchip:
#         ndates = mdates.date2num(data.time)
#         it = np.linspace(ndates[0], ndates[-1], num=ndates.size * 10)
#         # ft = mdates.num2date(it)
#         fx = pchip_interpolate(ndates, data.current_east.data, it)
#         fy = pchip_interpolate(ndates, data.current_north.data, it)
#     else:
#         it = mdates.date2num(data.time)
#         fx = data.current_east.data
#         fy = data.current_north.data
#     ftheta, frho = convert.cartesian_to_polar(fx, fy)
#     return it, ftheta, frho


# it, ftheta, frho = _get_surphys_currents_line(f_surphys, True)

# points = np.array([ftheta, frho]).T.reshape(-1, 1, 2)
# segments = np.concatenate([points[:-1], points[1:]], axis=1)

# fig = plt.figure(dpi=300)
# ax = fig.add_subplot(111, projection="polar")

# # Correct angle orientation
# ax.set_theta_zero_location("N")
# ax.set_theta_direction(-1)

# dc.plot.add_credit(ax)

# # ax.scatter(np.deg2rad(90), f_surphys.current_speed)
# ax.plot(ftheta, frho)

# norm = plt.Normalize(it.min(), it.max())
# lc = LineCollection(
#     segments, cmap="viridis", norm=norm, path_effects=[pe.Stroke(capstyle="round")]
# )
# # Set the values used for colormapping
# lc.set_array(it)
# lc.set_linewidth(5)
# line = ax.add_collection(lc)
