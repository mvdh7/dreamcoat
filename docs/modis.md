# MODIS satellite data

See the [NASA Ocean Color website](https://oceancolor.gsfc.nasa.gov/) for more information on the MODIS [Aqua](https://oceancolor.gsfc.nasa.gov/data/aqua/) and [Terra](https://oceancolor.gsfc.nasa.gov/data/terra/) instruments.  So far, only a small subset of the available data products can be accessed using dreamcoat, but this will be improved in future versions.  

## Set up credentials

You need to create a [NASA Earthdata account](https://urs.earthdata.nasa.gov/home) in order to download MODIS data.

Create a .netrc file and generate an AppKey as described [here](https://oceancolor.gsfc.nasa.gov/data/download_methods/) (check the *Download Methods* tab on that page).  Put the appkey in `~/.dreamcoat/nasa_appkey.dat`.

The download functions use `wget`, so they won't work natively on Windows.

## Download dataset

Dreamcoat downloads the daily L3-mapped observations from MODIS Aqua and Terra and combines them into a single dataset for each day, taking the mean of the two instruments.  To get a day or some consecutive days of data, use `get_days()`:

```python
import dreamcoat as dc

modis = dc.modis.get_days(
    date_min=None,
    date_max=None,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    appkey=None,
    delete_nc=False,
    filepath=".",
    resolution="9km",
)
```

  * Dates should be given in `%Y-%m-%d` format.  If left blank, yesterday is used for both.
  * Resolution can be `"4km"` or `"9km"`.
  * Files are downloaded as netCDFs into the `filepath`.

Only particulate inorganic carbon (PIC) data are downloaded.  These data are in mol/m<sup>3</sup>.

## Open daily datasets

To open a sequence of already-downloaded daily datasets and combine them into a single xarray Dataset, use `open_days()`:

```python
modis = dc.modis.open_days(
    date_min="1900-01-01",
    date_max="2300-01-01",
    filepath=".",
)
```
