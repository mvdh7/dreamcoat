import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import networkx as nx
from cartopy import crs as ccrs
from scipy import interpolate

from sys import path

epaths = ["/Users/matthew/github/neutralocean"]
for e in epaths:
    if e not in path:
        path.append(e)

from neutralocean import surface, label
import dreamcoat as dc

# %% Define and plot adjacency


class CruiseGraph(nx.Graph):
    def __init__(self, ctdz):
        super().__init__()
        self.ctdz = ctdz.copy()
        self._get_veronis_here()
        self._get_stations()
        self.levels_raw = {}
        self.levels = {}

    def _get_veronis_here(self):
        # Calculate Veronis densities with each station referenced to itself
        self.ctdz["veronis_here"] = np.nan
        for s, station in self.ctdz.groupby("station"):
            for i, row in station.iterrows():
                self.ctdz.loc[i, "veronis_here"] = label.veronis(
                    row.pressure,
                    station.salinity.values,
                    station.theta.values,
                    station.pressure.values,
                    eos="jmdfwg06",
                )
        self.ctdz["veronis_here"] -= 1000

    def _get_stations(self):
        self.stations = (
            self.ctdz[["station", "longitude", "latitude"]].groupby("station").mean()
        )
        self.stations["npts"] = (
            self.ctdz[["station", "longitude"]].groupby("station").count()
        )
        self.add_nodes_from(self.stations.index)

    def _get_edge_distances(self):
        distances = {}
        for _a, _b in self.edges:
            row_a = self.stations.loc[_a]
            row_b = self.stations.loc[_b]
            distances[(_a, _b)] = dc.maps.get_distance_geodesic(
                (row_a.longitude, row_a.latitude),
                (row_b.longitude, row_b.latitude),
            )
        nx.set_edge_attributes(self, distances, "distance")

    def add_edges_from(self, ebunch_to_add, **attr):
        assert np.all(
            np.isin(ebunch_to_add, self.stations.index)
        ), "All nodes must be in the existing list of stations!"
        super().add_edges_from(ebunch_to_add, **attr)
        self._get_edge_distances()
        self._get_grid()
        self._get_stp()

    def add_edge(self, u_of_edge, v_of_edge, **attr):
        assert (
            u_of_edge in self.stations.index and v_of_edge in self.stations.index
        ), "Both u_of_edge and v_of_edge must be in the existing list of stations!"
        super().add_edge(u_of_edge, v_of_edge, **attr)
        self._get_edge_distances()
        self._get_grid()
        self._get_stp()

    def _get_grid(self):
        a = []
        b = []
        distances = []
        for edge in self.edges:
            _a, _b = edge
            a.append((self.stations.index == _a).nonzero()[0][0])
            b.append((self.stations.index == _b).nonzero()[0][0])
            distances.append(self.edges[edge]["distance"])
        self.grid = {
            "edges": np.array([a, b]),
            "dist": np.array(distances),
            "distperp": np.ones_like(distances),
        }

    def _get_stp(self):
        stp_vars = ["salinity", "theta", "pressure"]
        stp = {}
        for v in stp_vars:
            stp[v] = np.full((self.stations.shape[0], self.stations.npts.max()), np.nan)
            for i, s in enumerate(self.stations.index):
                S = self.ctdz.station == s
                stp[v][i, : S.sum()] = self.ctdz[v][S].values
        self.stp = [stp[v] for v in stp_vars]

    def get_levels(self, station_ref):
        i_ref = (self.stations.index == station_ref).nonzero()[0][0]
        levels = np.full((self.stations.shape[0], self.stations.npts.max()), np.nan)
        for i, p_init in enumerate(
            self.ctdz[self.ctdz.station == station_ref].pressure
        ):
            levels[:, i] = surface.omega_surf(
                *self.stp, self.grid, pin_cast=i_ref, p_init=p_init, eos="jmdfwg06"
            )[2]
        self.levels_raw[station_ref] = levels
        self._smooth_and_interpolate(station_ref)
        self._get_ref_levels(station_ref)
        self._label_levels(station_ref)

    def _smooth_and_interpolate(self, station_ref):
        self.levels[station_ref] = np.full(self.levels_raw[station_ref].shape, np.nan)
        pressure_ref = self.get_pressure_ref(station_ref)
        for i in range(len(self.stations.index)):
            level = dc.plot.smooth_whittaker(self.levels_raw[station_ref][i])
            # Remove NaNs for interpolation and interpolate
            L = ~np.isnan(level)
            interp = interpolate.PchipInterpolator(
                pressure_ref[L], level[L], extrapolate=False
            )
            self.levels[station_ref][i] = interp(pressure_ref)

    def _get_ref_levels(self, station_ref):
        pressure_ref = self.get_pressure_ref(station_ref)
        self.ctdz["p_at_{}".format(station_ref)] = np.nan
        for i, s in enumerate(self.stations.index):
            L = ~np.isnan(self.levels[station_ref][i])
            interp = interpolate.PchipInterpolator(
                self.levels[station_ref][i][L], pressure_ref[L], extrapolate=False
            )
            S = self.ctdz.station == s
            self.ctdz.loc[S, "p_at_{}".format(station_ref)] = interp(
                self.ctdz.loc[S, "pressure"]
            )

    def _label_levels(self, station_ref):
        ref_veronis = self.ctdz[self.ctdz.station == station_ref][
            ["pressure", "veronis_here"]
        ]
        interp_veronis = interpolate.PchipInterpolator(
            ref_veronis.pressure, ref_veronis.veronis_here
        )
        self.ctdz["veronis_{}".format(station_ref)] = interp_veronis(
            self.ctdz["p_at_{}".format(station_ref)]
        )

    def get_pressure_ref(self, station_ref):
        i_ref = (graph.stations.index == station_ref).nonzero()[0][0]
        pressure_ref = graph.stp[2][i_ref]
        return pressure_ref

    def get_pos(self, crs):
        self.crs = crs
        self.pos = {
            i: self.crs.transform_points(
                ccrs.PlateCarree(), row.longitude, row.latitude
            )[0][:2]
            for i, row in self.stations.iterrows()
        }


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
graph = CruiseGraph(ctdz)
graph.add_edges_from(edges)


# Prepare edges for neutralocean surface calculation
extent = [-1, 7, 56.5, 64]
fcrs = ccrs.EquidistantConic(central_longitude=np.mean([extent[0], extent[1]]))
graph.get_pos(fcrs)

station_ref = 1
# station_ref = 86


graph.get_levels(station_ref)
pressure_ref = graph.get_pressure_ref(station_ref)


# %% Visualise surfaces at other stations
for i, s in enumerate(graph.stations.index):
    S = ctdz.station == s
    fig, ax = plt.subplots(dpi=300)
    ax.plot(pressure_ref, graph.levels_raw[station_ref][i], c="xkcd:black", lw=3)
    ax.plot(pressure_ref, graph.levels[station_ref][i], c="xkcd:strawberry", lw=3)
    ax.plot(
        "p_at_{}".format(station_ref),
        "pressure",
        data=graph.ctdz[S],
        c="xkcd:lime green",
        lw=0.8,
        zorder=3,
    )
    ax.scatter(pressure_ref, graph.levels_raw[station_ref][i], s=5, zorder=2)
    ax.set_title("Station {}".format(graph.stations.index[i]))
    ax.set_xlabel("Pressure at reference station / dbar")
    ax.set_ylabel("Connection depth at station / dbar")
    ax.set_aspect(1)
    fig.tight_layout()
    plt.show()
    plt.close()

# %% Visualise
fig, ax = plt.subplots(dpi=300)
ax.scatter(
    "veronis_here",
    "veronis_{}".format(station_ref),
    data=graph.ctdz,
    s=10,
    c="xkcd:sea blue",
    alpha=0.5,
    edgecolor="none",
)
ax.axline((27, 27), slope=1, c="k", lw=1)
ax.set_aspect(1)
ax.set_xlabel("Veronis anomaly per cast")
ax.set_ylabel("Veronis anomaly ref. to station {}".format(station_ref))

# %%
l = -1
l += 100
fig, ax = plt.subplots(dpi=300, subplot_kw={"projection": fcrs})
nx.draw_networkx(
    graph,
    pos=graph.pos,
    ax=ax,
    # Node marker properties
    node_size=0,
    linewidths=0,
    # Edge properties
    width=0.5,
    # Label properties
    font_size=8,
    with_labels=False,
)
sc = ax.scatter(
    "longitude",
    "latitude",
    data=graph.stations,
    transform=ccrs.PlateCarree(),
    s=20,
    c=graph.levels[station_ref][:, l],
    zorder=10,
    cmap="viridis_r",
    vmin=ctdz.loc[ctdz.station == station_ref].pressure.min(),
    vmax=ctdz.loc[ctdz.station == station_ref].pressure.max(),
)
plt.colorbar(sc, label="Neutral surface pressure / dbar")
ax.scatter(
    "longitude",
    "latitude",
    data=graph.stations.loc[station_ref],
    s=20,
    c="none",
    edgecolor="k",
    zorder=11,
    transform=ccrs.PlateCarree(),
)
ax.set_extent(extent, crs=ccrs.Geodetic())
ax.set_title(
    "Ref. pressure = {:.1f}".format(
        graph.ctdz.loc[graph.ctdz.station == station_ref].pressure.iloc[l]
    )
)
