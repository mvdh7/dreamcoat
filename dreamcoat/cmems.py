from datetime import date
import os


def download_cmems(
    date_min=None,
    date_max=None,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    filepath="",
    username=None,
    password=None,
):
    # if not username:
    #     if not keyring.get_password("dreamcoat.cmems", "dreamcoat.cmems.username"):
    #         keyring.set_password(
    #             "dreamcoat.cmems",
    #             "dreamcoat.cmems.username",
    #             input("Please enter your CMEMS username: "),
    #         )
    #     username = keyring.get_password("dreamcoat.cmems", "dreamcoat.cmems.username")
    # print(username)
    # if not keyring.get_password(
    #     "dreamcoat.cmems",
    #     keyring.get_password("dreamcoat.cmems", "dreamcoat.cmems.username"),
    # ):
    #     keyring.set_password(
    #         "dreamcoat.cmems",
    #         keyring.get_password("dreamcoat.cmems", "dreamcoat.cmems.username"),
    #         "\nPlease enter your CMEMS password: ",
    #     )
    date_min = "2022-10-24"
    date_max = "2022-11-02"
    if not date_min:
        date_min = date.today().strftime("%Y-%m-%d")
    if not date_max:
        date_max = date.today().strftime("%Y-%m-%d")
    os.system(
        "motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu "
        + "--service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS "
        + "--product-id global-analysis-forecast-phy-001-024 "
        + "--longitude-min {} ".format(longitude_min)
        + "--longitude-max {} ".format(longitude_max)
        + "--latitude-min {} ".format(latitude_min)
        + "--latitude-max {} ".format(latitude_max)
        + '--date-min "{} 12:00:00" '.format(date_min)
        + '--date-max "{} 12:00:00" '.format(date_max)
        + "--depth-min 0.494 --depth-max 0.4941 "
        + "--variable mlotst --variable so --variable thetao "
        + "--variable uo --variable vo --variable zos "
        + "--out-dir {} ".format(filepath)
        + "--out-name {}.nc ".format(filename)
        + "--user {} --pwd {}".format(username, password)
    )
