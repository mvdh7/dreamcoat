import pandas as pd
import dreamcoat as dc
from matplotlib import pyplot as plt

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
station_ref = 1
graph.get_surfaces(station_ref)

# Visualise surfaces at other stations
pressure_ref = graph.get_pressure_ref(station_ref)
for i, s in enumerate(graph.stations.index):
    S = ctdz.station == s
    fig, ax = plt.subplots(dpi=300)
    ax.plot(pressure_ref, graph.surfaces_raw[station_ref][i], c="xkcd:black", lw=3)
    ax.plot(pressure_ref, graph.surfaces[station_ref][i], c="xkcd:strawberry", lw=3)
    ax.plot(
        "p_at_{}".format(station_ref),
        "pressure",
        data=graph.ctdz[S],
        c="xkcd:lime green",
        lw=0.8,
        zorder=3,
    )
    ax.scatter(pressure_ref, graph.surfaces_raw[station_ref][i], s=5, zorder=2)
    ax.set_title("Station {}".format(graph.stations.index[i]))
    ax.set_xlabel("Pressure at reference station / dbar")
    ax.set_ylabel("Connection depth at station / dbar")
    ax.set_aspect(1)
    fig.tight_layout()
    plt.show()
    plt.close()
