# MODIS satellite data

See the [NASA Ocean Color website](https://oceancolor.gsfc.nasa.gov/) for more information on the MODIS [Aqua](https://oceancolor.gsfc.nasa.gov/data/aqua/) and [Terra](https://oceancolor.gsfc.nasa.gov/data/terra/) instruments.  So far, only a small subset of the available data products can be accessed using dreamcoat, but this will be improved in future versions.  

## Download dataset

Dreamcoat downloads the daily datasets from MODIS Aqua and Terra and combines them into a single dataset for each day, taking the mean of the two instruments.  To get a day or some consecutive days of data, use `get_modis_days()`:

```python
import dreamcoat as dc

modis = dc.modis.get_modis_days(
    date_min=None,
    date_max=None,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    resolution="9km",
)
```

  * Dates should be given in `%Y-%m-%d` format.  If left blank, yesterday is used for both.
  * Resolution can be `"4km"` or `"9km"`.

Only particulate inorganic carbon (PIC) data are downloaded.  These data are in mol/m<sup>3</sup>.
