version_number = "0.1"
version_colour = "red"
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
