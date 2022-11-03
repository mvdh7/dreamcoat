import numpy as np
import dreamcoat as dc


# boxros-6qacgu-bUqmyr

lims = dict(
    date_min="2021-05-20",
    date_max="2021-06-08",
    longitude_min=-10,
    longitude_max=10,
    latitude_min=55,
    latitude_max=65,
)
# dc.cmems.download_surphys(
# filepath="tests/data/",
# **lims,
# )
phys = dc.cmems.open_surphys(
    filepath="tests/data/",
    **lims,
)


speed_knots = 9
timestep_hours = 24
ship_distance = (
    dc.convert.knots_to_kph(speed_knots) * timestep_hours * np.array([1, 2, 3])
)

#%%
for i in range(20):
    fig, ax = dc.plot.surphys(
        # phys.mean('time'),
        phys.isel(time=i),
        "current_speed",
        land_visible=False,
        # ship_longitude=0,
        # ship_latitude=-35,
        # ship_distance=ship_distance,
        # quiver_coarsen=20,
        # quiver_alpha=0.1,
        map_extent=[-3, 10, 57, 64],
        # vmin=0,
        # vmax=1,
    )
