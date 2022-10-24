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


phys = dc.cmems.open_cmems(
    filepath="tests/data/",
    longitude_min=-70,
    longitude_max=30,
    latitude_min=-50,
    latitude_max=-20,
)
