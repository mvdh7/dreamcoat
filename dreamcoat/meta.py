import os

version_number = "0.3"
version_colour = "green"
__version__ = "{}-{}".format(version_number, version_colour)


def hello():
    print(
        r"""
    ____    ____    U _____ u    _      __  __     ____   U  ___ u    _       _____   
   |  _"\U |  _"\ u \| ___"|/U  /"\  uU|' \/ '|uU /"___|   \/"_ \/U  /"\  u  |_ " _|  
  /| | | |\| |_) |/  |  _|"   \/ _ \/ \| |\/| |/\| | u     | | | | \/ _ \/     | |    
  U| |_| |\|  _ <    | |___   / ___ \  | |  | |  | |/__.-,_| |_| | / ___ \    /| |\   
   |____/ u|_| \_\   |_____| /_/   \_\ |_|  |_|   \____|\_)-\___/ /_/   \_\  u |_|U   
    |||_   //   \\_  <<   >>  \\    >><<,-,,-.   _// \\      \\    \\    >>  _// \\_  
   (__)_) (__)  (__)(__) (__)(__)  (__)(./  \.) (__)(__)    (__)  (__)  (__)(__) (__)
   """
    )
    print(
        "{:^88}".format("The joyful oceanographic seagoing expedition planning helper")
    )
    print(
        "{:^88}".format(
            "{} version ({})".format(version_colour.title(), version_number)
        )
    )


def get_dat_data(filename):
    """Find data from a file in the .dreamcoat folder, which should be in the user's
    home path, as returned by `os.path.expanduser('~')`.  If the file is not found, then
    the user is prompted to enter a value instead.

    Parameters
    ----------
    filename : str
        The name of the file, excluding the .dat extension.

    Returns
    -------
    str
        The value stored in the .dat file or provided by the user.
    """
    fn = os.sep.join((os.path.expanduser("~"), ".dreamcoat", "{}.dat".format(filename)))
    try:
        with open(fn, "r") as f:
            data = f.read().splitlines()[0]
    except FileNotFoundError:
        print("File {} not found!".format(fn))
        data = input("Please enter a value for {}: ".format(filename))
    return data
