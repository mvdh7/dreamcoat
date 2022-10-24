from getpass import getpass
import keyring
import os


def download_cmems(
    filepath="",
    username=None,
    password=None,
):
    if not username:
        if not keyring.get_password("dreamcoat.cmems", "dreamcoat.cmems.username"):
            keyring.set_password(
                "dreamcoat.cmems",
                "dreamcoat.cmems.username",
                input("Please enter your CMEMS username: "),
            )
        username = keyring.get_password("dreamcoat.cmems", "dreamcoat.cmems.username")
    print(username)
    # if not keyring.get_password(
    #     "dreamcoat.cmems",
    #     keyring.get_password("dreamcoat.cmems", "dreamcoat.cmems.username"),
    # ):
    #     keyring.set_password(
    #         "dreamcoat.cmems",
    #         keyring.get_password("dreamcoat.cmems", "dreamcoat.cmems.username"),
    #         "\nPlease enter your CMEMS password: ",
    #     )
    # os.system(
    #     "motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu "
    #     + "--service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS "
    #     + "--product-id global-analysis-forecast-phy-001-024 "
    #     + "--longitude-min -70 "
    #     + "--longitude-max 30 "
    #     + "--latitude-min -60 "
    #     + "--latitude-max -20 "
    #     + '--date-min "2022-10-24 12:00:00" '
    #     + '--date-max "2022-11-02 12:00:00" '
    #     + "--depth-min 0.494 --depth-max 0.4941 "
    #     + "--variable mlotst --variable so --variable thetao "
    #     + "--variable uo --variable vo --variable zos "
    #     + "--out-dir {} ".format(filepath)
    #     + "--out-name test.nc "
    #     + "--user {} --pwd {}".format(username, password)
    # )
