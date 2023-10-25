# Goal: distance-to-shore calculation for NWESDAP

import dreamcoat as dc
import numpy as np
from matplotlib import pyplot as plt
from cartopy import crs as ccrs, feature as cfeature

# Build VP tree
vptree_lats = np.array([42, 60])
vptree_lons = np.array([-12, 21])
vpt = dc.build_vptree(vptree_lats, vptree_lons, resolution="10m")

# Use VP tree
vlat = np.linspace(10, 80, 70)
vlon = np.linspace(-40, 60, 100)
lon, lat = np.meshgrid(vlon, vlat)
distance = np.full_like(lon, np.nan)

for a in range(len(vlat)):
    print(a)
    for o in range(len(vlon)):
        distance[a, o] = vpt.get_nearest_neighbor((vlon[o], vlat[a]))[0]

# Plot result
fig, ax = plt.subplots(dpi=300, subplot_kw={"projection": ccrs.PlateCarree()})
sc = ax.scatter(
    lon.ravel(),
    lat.ravel(),
    c=distance.ravel(),
    s=12,
    marker="s",
    edgecolor="none",
    transform=ccrs.PlateCarree(),
    cmap="turbo",
)
ax.add_feature(cfeature.COASTLINE)
ax.plot(
    vptree_lons[[0, 0, 1, 1, 0]],
    vptree_lats[[0, 1, 1, 0, 0]],
    transform=ccrs.PlateCarree(),
    c="magenta",
)
plt.colorbar(sc)
