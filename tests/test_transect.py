import pandas as pd
import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
import dreamcoat as dc
from neutralocean.traj import neutral_trajectory
from neutralocean.surface import omega_surf
from neutralocean.label import veronis

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
transect = transects[1]
cvar = "dic"

# kwargs
extra_fraction = 0.05

# Subset to transect
tdeep = deep[deep.station.isin(transect)].copy()
tctdz = ctdz[ctdz.station.isin(transect)].copy()
tstations = (
    tdeep[["station", "latitude", "longitude"]]
    .groupby("station")
    .mean()
    .reindex(transect)
)
tstations["distance"] = dc.maps.get_route_distance(
    tstations[["longitude", "latitude"]].values.T
)
tstations["npts_ctdz"] = tctdz[["station", "longitude"]].groupby("station").count()
tdeep["transect_distance"] = tstations.loc[tdeep.station].distance.values
tctdz["transect_distance"] = tstations.loc[tctdz.station].distance.values

# Make and extend route
route_lon, route_lat, route_distance = dc.maps.extend_route(
    tstations.longitude.values, tstations.latitude.values, extra_fraction=extra_fraction
)

# Calculate Veronis densities with each station referenced to itself
tctdz["veronis_here"] = np.nan
for s, station in tctdz.groupby("station"):
    for i, row in station.iterrows():
        tctdz.loc[i, "veronis_here"] = veronis(
            row.pressure,
            station.salinity.values,
            station.theta.values,
            station.pressure.values,
            eos="jmdfwg06",
        )
tctdz["veronis_here"] -= 1000


def get_stp(tctdz, tstations):
    stp_vars = ["salinity", "theta", "pressure"]
    stp = {}
    for v in stp_vars:
        stp[v] = np.full((tstations.shape[0], tstations.npts_ctdz.max()), np.nan)
        for i, s in enumerate(tstations.index):
            S = tctdz.station == s
            stp[v][i, : S.sum()] = tctdz[v][S].values
    stp = np.array([stp[v] for v in stp_vars])
    return stp


s_iref = 0  # the iloc index of the reference station

# Get neutral trajectories
stp = get_stp(tctdz, tstations)
nts = []
nsurf = []
for p0 in tctdz[tctdz.station == tstations.index[s_iref]].pressure:
    nts.append(neutral_trajectory(*stp, p0, eos="jmdfwg06")[2])
    # nsurf.append(omega_surf(
    #     *stp, grid, pin_cast=i_ref, p_init=p_init, eos="jmdfwg06"
    # )[2])
nts = np.array(nts).T

# TODO I think that using omega surfaces will be better than neutral trajectories here
# because we can assign distances and more easily switch between different stations
# as the reference (with neutral trajectories, the first station is always the ref.)
# BUT to implement this it would be much more efficient to first tidy up and then use
# the graph tools currently in surfaces.py --- either make a new package or just throw
# them into dreamcoat for now (latter probably easier).

# Smooth the neutral trajectories so they (hopefully) only increase with pressure
# (this step does not currently guarantee monotonicity)
window = 7
for i, s in enumerate(tstations.index):
    S = tctdz.station == tstations.index[s_iref]
    tctdz.loc[S, "p_at_{}".format(s)] = nts[i]
    tctdz.loc[S, "p_at_{}_s".format(s)] = (
        tctdz.loc[S, "p_at_{}".format(s)]
        .rolling(window, min_periods=5, win_type="gaussian")
        .mean(std=1)  # change the value of std to adjust the degree of smoothing
    ).shift(-int(np.floor(window / 2)))

# Invert trajectories (get pressure at reference station corresponding to others)
tctdz["p_from_{}".format(tstations.index[s_iref])] = np.nan
for s in tstations.index:
    L = (tctdz.station == tstations.index[s_iref]) & tctdz[
        "p_at_{}_s".format(s)
    ].notnull()
    interp = interpolate.PchipInterpolator(
        tctdz[L].pressure, tctdz[L]["p_at_{}_s".format(s)], extrapolate=False
    )
    S = tctdz.station == s
    tctdz.loc[S, "p_from_{}".format(tstations.index[s_iref])] = interp(
        tctdz[S].pressure
    )

# Visualise trajectory smoothing
fig, ax = plt.subplots(dpi=300)
for s in tstations.index:
    ax.scatter(
        "p_at_{}".format(tstations.index[s_iref]),
        "p_at_{}".format(s),
        data=tctdz,
        label=s,
        s=10,
        alpha=0.2,
    )
    ax.plot(
        "p_at_{}_s".format(tstations.index[s_iref]),
        "p_at_{}_s".format(s),
        data=tctdz,
        label="",
    )
ax.set_xlabel("Pressure at station {} / dbar".format(tstations.index[s_iref]))
ax.set_ylabel("Neutral trajectory pressure / dbar")
ax.legend()
ax.grid(alpha=0.2)

# Visualise inversion
fig, ax = plt.subplots(dpi=300)
for s in tstations.index:
    ax.plot(
        "p_at_{}_s".format(tstations.index[s_iref]),
        "p_at_{}_s".format(s),
        data=tctdz,
        label="",
    )
    ax.plot(
        "pressure",
        "p_from_{}".format(tstations.index[s_iref]),
        data=tctdz[tctdz.station == s],
        label="",
        ls="--",
    )
ax.set_xlabel("Pressure at station {} / dbar".format(tstations.index[s_iref]))
ax.set_ylabel("Neutral trajectory pressure / dbar")
ax.grid(alpha=0.2)

# %% Visualise transect
fig, ax = plt.subplots(dpi=300)
# ax.scatter("transect_distance", "pressure", data=tdeep, c=cvar)
for s, station in tctdz.groupby("station"):
    ax.scatter(
        "transect_distance",
        # "pressure",
        # c="p_from_{}".format(tstations.index[s_iref]),
        "p_from_{}".format(tstations.index[s_iref]),
        c="theta",
        data=station,
    )
ax.invert_yaxis()
ax.set_xlim(route_distance[0], route_distance[-1])
ax.set_xlabel("Transect distance / km")
ax.set_ylabel("Pressure / dbar")
