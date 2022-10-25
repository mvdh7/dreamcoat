import numpy as np
import dreamcoat as dc


# dc.cmems.download_cmems(
#     filepath="tests/data/",
#     longitude_min=-70,
#     longitude_max=30,
#     latitude_min=-60,
#     latitude_max=-20,
#     date_min="2022-10-24",
#     date_max="2022-11-02",
# )

# boxros-6qacgu-bUqmyr

lat_lon_lims = dict(
    longitude_min=-70,
    longitude_max=30,
    latitude_min=-40,
    latitude_max=-20,
)
phys = dc.cmems.open_surphys(
    filepath="tests/data/",
    **lat_lon_lims,
    date_min="2022-10-24",
    date_max="2022-10-24",
)


speed_knots = 9
timestep_hours = 24
ship_distance = (
    dc.convert.knots_to_kph(speed_knots) * timestep_hours * np.array([1, 2, 3])
)


fig, ax = dc.plot.surphys(
    phys,
    "mld",
    land_visible=False,
    ship_longitude=0,
    ship_latitude=-35,
    ship_distance=ship_distance,
    quiver_coarsen=20,
    quiver_alpha=0.1,
    map_extent=[-5, 5, -39, -31],
    # vmin=34,
    # vmax=37,
)
