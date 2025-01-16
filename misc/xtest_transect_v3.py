import pandas as pd
import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
import dreamcoat as dc

# args
deep = pd.read_parquet("tests/data/deep.parquet")
ctdz = pd.read_parquet("tests/data/ctdz.parquet")
transects = {}
transects[0] = [1, 59, 112, 81, 86]
transects[1] = [8, 46, 63, 93]
transects[2] = [1, 8, 19, 24, 27]
transects[3] = [35, 40, 46, 50, 59]
transects[4] = [70, 63, 74, 112, 117]
transects[5] = [93, 89, 86, 98]
transect_stations = transects[3]

# kwargs
extra_fraction = 0.05

# Subset to transect
tdeep = deep[deep.station.isin(transect_stations)].copy()
tctdz = ctdz[ctdz.station.isin(transect_stations)].copy()
tstations = (
    tdeep[["station", "latitude", "longitude"]]
    .groupby("station")
    .mean()
    .reindex(transect_stations)
)
tstations["distance"] = dc.maps.get_route_distance(
    tstations[["longitude", "latitude"]].values.T
)
tdeep["transect_distance"] = tstations.loc[tdeep.station].distance.values
tctdz["transect_distance"] = tstations.loc[tctdz.station].distance.values

# Make and extend route
route_lon, route_lat, route_distance = dc.maps.extend_route(
    tstations.longitude.values, tstations.latitude.values, extra_fraction=extra_fraction
)

# Get neutral trajectories
transect = dc.CruiseGraph(tctdz)
for i in range(1, len(transect_stations)):
    transect.add_edge(transect_stations[i - 1], transect_stations[i])

from neutralocean.traj import neutral_trajectory

# %%
tc = transect.ctdz  # just for convenience
s = 0
S = tc.station == transect_stations[s]
tc.loc[S, "veronis"] = tc[S].veronis_here.copy()
for s in range(len(transect_stations) - 1):
    S = tc.station == transect_stations[s]
    traj = []
    for p in tc[S].pressure:
        traj.append(
            neutral_trajectory(
                transect.stp[0][s : s + 2],
                transect.stp[1][s : s + 2],
                transect.stp[2][s : s + 2],
                p,
                eos="jmdfwg06",
            )[2][1]
        )
    traj = np.array(traj)
    traj = dc.plot.smooth_whittaker(traj)
    L = ~np.isnan(traj)
    interp_traj = interpolate.PchipInterpolator(
        tc[S].pressure[L], traj[L], extrapolate=False
    )
    traj = interp_traj(tc[S].pressure.values)
    # Invert trajectory to get Veronis labels from previous station
    L = ~np.isnan(traj)
    interp_invert = interpolate.PchipInterpolator(
        traj[L], tc[S].pressure[L], extrapolate=False
    )
    T = tc.station == transect_stations[s + 1]
    tc.loc[T, "p_at_{}".format(transect_stations[s])] = interp_invert(tc[T].pressure)
    interp_veronis = interpolate.PchipInterpolator(
        tc[S].pressure, tc[S].veronis, extrapolate=False
    )
    tc.loc[T, "veronis_{}".format(transect_stations[s])] = interp_veronis(
        tc[T]["p_at_{}".format(transect_stations[s])]
    )
    tc.loc[T, "veronis"] = tc[T]["veronis_{}".format(transect_stations[s])].copy()
    # Find where the next station along the transect doesn't have Veronis labels
    # assigned from the previous station (at the top and bottom)
    veronis_ix_range = tc[T].veronis.notnull().values.nonzero()[0][[0, -1]]
    # ^ this contains the [first, last] iloc with Veronis values
    # Calculate the difference between the in situ Veronis label and the labels
    # assigned from the previous station at the top and bottom of the assignment range
    veronis_offset = (tc.veronis - tc.veronis_here)[T].iloc[veronis_ix_range].values
    # Fill in the unassigned Veronis labels with the in situ values adjusted by the
    # offset, so that the Veronis curve joins smoothly and remains monotonic
    veronis_next = tc[T].veronis.values
    veronis_here = tc[T].veronis_here.values
    veronis_next[: veronis_ix_range[0]] = (
        veronis_here[: veronis_ix_range[0]] + veronis_offset[0]
    )
    veronis_next[veronis_ix_range[1] :] = (
        veronis_here[veronis_ix_range[1] :] + veronis_offset[1]
    )
    tc.loc[T, "veronis"] = veronis_next

# %% Assign Veronis labels to deep
for s in transect_stations:
    Stc = tc.station == s
    interp_veronis = interpolate.PchipInterpolator(
        tc[Stc].pressure, tc[Stc].veronis, extrapolate=True
    )
    S = tdeep.station == s
    tdeep.loc[S, "veronis"] = interp_veronis(tdeep[S].pressure)
    # Set any extrapolated values that are beyond range to min/max
    Sx = S & (tdeep[S].veronis < tc[Stc].veronis.min())
    if Sx.any():
        tdeep.loc[Sx, "veronis"] = tc[Stc].veronis.min()
    Sx = S & (tdeep[S].veronis > tc[Stc].veronis.max())
    if Sx.any():
        tdeep.loc[Sx, "veronis"] = tc[Stc].veronis.max()

# %% Visualise transect
cvar = "salinity"
# cvar = "p_at_{}".format(transect_stations[i_ref])
fig, ax = plt.subplots(dpi=300)
# ax.scatter("transect_distance", "pressure", data=tdeep, c=cvar)
for s, station in transect.ctdz.groupby("station"):
    fsc = ax.scatter(
        "transect_distance",
        "veronis",
        c=cvar,
        cmap="viridis_r",
        data=station,
        vmin=transect.ctdz[cvar].min(),
        vmax=transect.ctdz[cvar].max(),
    )
plt.colorbar(fsc, label=cvar)
ax.invert_yaxis()
ax.set_xlim(route_distance[0], route_distance[-1])
ax.set_xlabel("Transect distance / km")
ax.set_ylabel("Veronis anomaly")

# %% Generate interpolation
skip_station = transect_stations[2]
kernel = dict(
    kernel="multiquadric",
    epsilon=100,
)

# Get interpolator to find cvar on a distance-Veronis grid

# # With the high-res data
# L = np.isin(transect.ctdz.station, transect_stations)
# x = transect.ctdz[L].transect_distance.values
# y = transect.ctdz[L].veronis.values
# z = transect.ctdz[L][cvar].values

# Or with low-res data
x = tdeep.transect_distance.values
y = tdeep.veronis.values
z = tdeep[cvar].values

# Either way, carry on here
L = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(z)
L &= tdeep.station != skip_station  # to see effect of missing data
x = x[L]
y = y[L]
z = z[L]
xscale_cvar = 0.005
xy = np.array([x * xscale_cvar, y]).T
interp_cvar = interpolate.RBFInterpolator(xy, z, **kernel)
# interp_cvar = interpolate.LinearNDInterpolator(xy, z)

# Visualise cvar interpolator
gx = np.linspace(route_distance[0], route_distance[-1], num=500)
gy = np.linspace(np.min(y), np.max(y), num=500)
gx, gy = np.meshgrid(gx, gy)
gz = interp_cvar(np.array([gx.ravel() * xscale_cvar, gy.ravel()]).T)
gz = np.reshape(gz, gx.shape)
fig, ax = plt.subplots(dpi=300)
fc = ax.contourf(
    gx,
    gy,
    gz,
    80,
    vmin=transect.ctdz[cvar].min(),
    vmax=transect.ctdz[cvar].max(),
)
fs = ax.scatter(
    x,
    y,
    c=z,
    s=3,
    vmin=transect.ctdz[cvar].min(),
    vmax=transect.ctdz[cvar].max(),
)
plt.colorbar(fs, label=cvar)
ax.set_title("Interpolate cvar on distance and Veronis")
ax.set_xlabel("Transect distance / km")
ax.set_ylabel("Veronis anomaly")

# Get interpolator to find Veronis on a distance-pressure grid
L = np.isin(transect.ctdz.station, transect_stations)
L &= transect.ctdz.station != skip_station  # to see effect of missing data
x_v = transect.ctdz[L].transect_distance.values
y_v = transect.ctdz[L].pressure.values
z_v = transect.ctdz[L].veronis.values
L = ~np.isnan(x_v) & ~np.isnan(y_v) & ~np.isnan(z_v)
x_v = x_v[L]
y_v = y_v[L]
z_v = z_v[L]
xscale_veronis = 0.5
xy_v = np.array([x_v * xscale_veronis, y_v]).T
interp_veronis = interpolate.RBFInterpolator(xy_v, z_v, **kernel)
# interp_veronis = interpolate.LinearNDInterpolator(xy_v, z_v)

# Create grid and predict cvar on it
gx_v = np.linspace(route_distance[0], route_distance[-1], num=500)
gy_v = np.linspace(0, np.max(y_v) + 10, num=500)
gx_v, gy_v = np.meshgrid(gx_v, gy_v)
gz_v = interp_veronis(np.array([gx_v.ravel() * xscale_veronis, gy_v.ravel()]).T)
gz_v = np.reshape(gz_v, gx_v.shape)

# Visualise Veronis interpolation
fig, ax = plt.subplots(dpi=300)
fc = ax.contourf(
    gx_v,
    gy_v,
    gz_v,
    256,
    vmin=transect.ctdz.veronis.min(),
    vmax=transect.ctdz.veronis.max(),
)
ax.contour(
    gx_v,
    gy_v,
    gz_v,
    [24, 25, 26, *np.arange(27, 28.5, 0.1)],
    colors="k",
    linewidths=0.8,
    alpha=0.1,
    vmin=transect.ctdz.veronis.min(),
    vmax=transect.ctdz.veronis.max(),
)
fs = ax.scatter(
    x_v,
    y_v,
    c=z_v,
    s=3,
    vmin=transect.ctdz.veronis.min(),
    vmax=transect.ctdz.veronis.max(),
)
plt.colorbar(fs, label="Veronis")
ax.set_ylim(y_v.max(), 0)
ax.set_title("Interpolate Veronis on distance and pressure")
ax.set_xlabel("Transect distance / km")
ax.set_ylabel("Pressure / dbar")

# Combine interpolators to get final cvar field
# Predict cvar on Veronis grid
gc = interp_cvar(np.array([gx_v.ravel() * xscale_cvar, gz_v.ravel()]).T)
gc = np.reshape(gc, gx_v.shape)

# Visualise Veronis interpolation
fig, ax = plt.subplots(dpi=300)
fc = ax.contourf(
    gx_v,
    gy_v,
    gc,
    256,
    vmin=transect.ctdz[cvar].min(),
    vmax=transect.ctdz[cvar].max(),
)
ax.contour(
    gx_v,
    gy_v,
    gz_v,
    [24, 25, 26, *np.arange(27, 28.5, 0.1)],
    colors="xkcd:dark",
    linewidths=0.8,
    alpha=0.3,
)
fs = ax.scatter(
    "transect_distance",
    "pressure",
    c=cvar,
    data=transect.ctdz,
    s=3,
    vmin=transect.ctdz[cvar].min(),
    vmax=transect.ctdz[cvar].max(),
)
plt.colorbar(fs, label=cvar)
ax.set_ylim(gy_v.max(), 0)
ax.set_title("Interpolate cvar on distance and pressure via Veronis")
ax.set_xlabel("Transect distance / km")
ax.set_ylabel("Pressure / dbar")


# %% Interpolate cvar directly on pressure

# # With high-res data
# L = np.isin(transect.ctdz.station, transect_stations)
# x_d = transect.ctdz[L].transect_distance.values
# y_d = transect.ctdz[L].pressure.values
# z_d = transect.ctdz[L][cvar].values

# Or with low-res data
x_d = tdeep.transect_distance.values
y_d = tdeep.pressure.values
z_d = tdeep[cvar].values

L = ~np.isnan(x_d) & ~np.isnan(y_d) & ~np.isnan(z_d)
L &= tdeep.station != skip_station  # to see effect of missing data
x_d = x_d[L]
y_d = y_d[L]
z_d = z_d[L]
xscale_direct = 0.5
xy_d = np.array([x_d * xscale_direct, y_d]).T
interp_direct = interpolate.RBFInterpolator(xy_d, z_d, **kernel)
# interp_direct = interpolate.LinearNDInterpolator(xy_d, z_d)

# Create grid and predict cvar on it
gx_d = np.linspace(route_distance[0], route_distance[-1], num=500)
gy_d = np.linspace(0, np.max(y_d) + 10, num=500)
gx_d, gy_d = np.meshgrid(gx_d, gy_d)
gz_d = interp_direct(np.array([gx_d.ravel() * xscale_direct, gy_d.ravel()]).T)
gz_d = np.reshape(gz_d, gx_d.shape)

# Visualise interpolation
fig, ax = plt.subplots(dpi=300)
fc = ax.contourf(
    gx_d,
    gy_d,
    gz_d,
    256,
    vmin=transect.ctdz[cvar].min(),
    vmax=transect.ctdz[cvar].max(),
)
ax.contour(
    gx_v,
    gy_v,
    gz_v,
    [24, 25, 26, *np.arange(27, 28.5, 0.1)],
    colors="xkcd:dark",
    linewidths=0.8,
    alpha=0.3,
    vmin=transect.ctdz.veronis.min(),
    vmax=transect.ctdz.veronis.max(),
)
fs = ax.scatter(
    # x_d,
    # y_d,
    # c=z_d,
    transect.ctdz.transect_distance,
    transect.ctdz.pressure,
    c=transect.ctdz[cvar],
    s=3,
    vmin=transect.ctdz[cvar].min(),
    vmax=transect.ctdz[cvar].max(),
)
plt.colorbar(fs, label=cvar)
ax.set_ylim(y_d.max(), 0)
ax.set_title("Interpolate cvar on distance and pressure directly")
ax.set_xlabel("Transect distance / km")
ax.set_ylabel("Pressure / dbar")

# TODO make some way to quantify differences between in situ samples and interpolated
# fields, for meaningful testing and comparisons.
# Looks like via Veronis does work better with missing stations!
