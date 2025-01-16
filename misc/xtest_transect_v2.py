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
transect_stations = transects[1]
cvar = "dic"

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

# Get omega surfaces
transect = dc.CruiseGraph(tctdz)
for i in range(1, len(transect_stations)):
    transect.add_edge(transect_stations[i - 1], transect_stations[i])
transect.get_surfaces_all()

# %% Visualise transect
i_ref = -1
i_ref += 1
# cvar = "theta"
cvar = "p_at_{}".format(transect_stations[i_ref])
fig, ax = plt.subplots(dpi=300)
# ax.scatter("transect_distance", "pressure", data=tdeep, c=cvar)
for s, station in transect.ctdz.groupby("station"):
    fsc = ax.scatter(
        "transect_distance",
        "p_at_{}".format(transect_stations[i_ref]),
        # "pressure",
        c=cvar,
        cmap="viridis_r",
        data=station,
        vmin=transect.ctdz[cvar].min(),
        vmax=transect.ctdz[cvar].max(),
    )
plt.colorbar(fsc)
ax.invert_yaxis()
ax.set_xlim(route_distance[0], route_distance[-1])
ax.set_xlabel("Transect distance / km")
ax.set_ylabel("Pressure / dbar")

# %% Generate interpolation
i_ref = 3
# cvar = "p_at_{}".format(transect_stations[i_ref])
cvar = "veronis_{}".format(transect_stations[i_ref])

L = np.isin(transect.ctdz["station"], transect_stations)
x = transect.ctdz[L]["transect_distance"].values
y = transect.ctdz[L]["pressure"].values
z = transect.ctdz[L][cvar].values
L = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(z)
x = x[L]
y = y[L]
z = z[L]
xscale = 0.1
xy = np.array([x * xscale, y]).T
# interp = interpolate.RBFInterpolator(xy, z, kernel="linear")
interp = interpolate.LinearNDInterpolator(xy, z)

# Create grid and predict cvar on it
gx = np.linspace(route_distance[0], route_distance[-1], num=200)
gy = np.linspace(0, np.max(y) + 10, num=200)
gx, gy = np.meshgrid(gx, gy)
gz = interp(np.array([gx.ravel() * xscale, gy.ravel()]).T)
gz = np.reshape(gz, gx.shape)


fig, ax = plt.subplots(dpi=300)
fc = ax.contourf(
    gx,
    gy,
    gz,
    80,
    vmin=transect.ctdz.veronis_here.min(),
    vmax=transect.ctdz.veronis_here.max(),
    # vmin=transect.ctdz[cvar].min(),
    # vmax=transect.ctdz[cvar].max(),
    # vmin=0,
    # vmax=transect.ctdz.pressure.max(),
)
fs = ax.scatter(
    x,
    y,
    c=z,
    vmin=transect.ctdz.veronis_here.min(),
    vmax=transect.ctdz.veronis_here.max(),
    s=3,
    # vmin=transect.ctdz[cvar].min(),
    # vmax=transect.ctdz[cvar].max(),
    # vmin=0,
    # vmax=transect.ctdz.pressure.max(),
)
plt.colorbar(fs)
ax.set_ylim(transect.ctdz.pressure.max() + 10, 0)

# %% Make final Veronis labels
# transects[1] = [8, 46, 63, 93]
tc = transect.ctdz  # just for convenience

# Start with the Veronis labels propagated from the first station on the transect
s = 0
tc["veronis"] = tc["veronis_{}".format(transect_stations[s])].copy()

# %% Find where the next station along the transect doesn't have Veronis labels assigned
# from the previous station (at the top and bottom)
s += 1
S = tc.station == transect_stations[s]
veronis_ix_range = tc[S].veronis.notnull().values.nonzero()[0][[0, -1]]
# ^ this contains the [first, last] iloc with Veronis values

# Calculate the difference between the in situ Veronis label and the labels assigned
# from the previous station at the top and bottom of the assignment range
veronis_offset = (tc.veronis - tc.veronis_here)[S].iloc[veronis_ix_range].values

# Fill in the unassigned Veronis labels with the in situ values adjusted by the offset,
# so that the Veronis curve joins smoothly and remains monotonic
veronis_next = tc[S].veronis.values
veronis_here = tc[S].veronis_here.values
veronis_next[: veronis_ix_range[0]] = (
    veronis_here[: veronis_ix_range[0]] + veronis_offset[0]
)
veronis_next[veronis_ix_range[1] :] = (
    veronis_here[veronis_ix_range[1] :] + veronis_offset[1]
)
tc.loc[S, "veronis"] = veronis_next

# # Reassign Veronis labels for the next station along (third station)
# interp = interpolate.PchipInterpolator(tc[S].pressure, tc[S].veronis)
# S = tc.station == transect_stations[s + 1]
# tc.loc[S, "veronis"] = interp(tc[S]["p_at_{}".format(transect_stations[1])])

# %% Visualise
fs = transect_stations[s]
fdata = tc[tc.station == fs]
fig, ax = plt.subplots()
ax.plot("pressure", "veronis_here", data=fdata)
ax.plot("pressure", "veronis", data=fdata)
ax.plot("pressure", "veronis_{}".format(transect_stations[s - 1]), data=fdata)
ax.plot("pressure", "veronis_8", data=fdata)
