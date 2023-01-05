# Working with CMEMS datasets

You can use dreamcoat to download, process and plot daily surface ocean reanalysis / forecast data from the CMEMS datasets

  * `GLOBAL_ANALYSIS_FORECAST_PHY_001_024` and 
  * `GLOBAL_ANALYSIS_FORECAST_BIO_001_028`.

These two datasets are collected for the same space and time ranges and merged into a single xarray Dataset for further processing and plotting.  The biogeochemical dataset comes at a lower resolution than the physical one so the former is nearest-neighbour interpolated onto the latter's latitude-longitude grid.

New files are only downloaded if they don't already exist locally and only the specified space and time ranges are downloaded.

## Download or open dataset

```python
import dreamcoat as dc

surface = dc.cmems.open_surface(
    filepath="",
    date_min=None,
    date_max=None,
    latitude_min=-90,
    latitude_max=90,
    longitude_min=-180,
    longitude_max=180,
    username=None,
    password=None,
)
```

  * The `filepath` specifies where downloaded files should be saved.
  * Dates should be given in `%Y-%m-%d` format.  If left blank, today is used for both.
  * Username and password should be your login details for the CMEMS data portal.  **Don't commit these to a public repo!**

## Plotting functions

### Surface maps

#### At a single time point

`surface_map` draws a map of a single time point from the CMEMS dataset.  You need to extract the single time slice yourself before running the plotting function, as shown below.

```python
surface_time_slice = surface.isel(0)  # for example
surface_time_slice = surface.mean('time')  # alternative example

fig, ax = dc.plot.surface_map(
    surface_time_slice,
    fvar,  # string specifying which variable to plot
    ax=None,
    color_zoom_factor=0.9,
    dpi=150,
    figsize=[6.4, 4.8],
    land_visible=True,
    longitude_fmt=":03.0f",
    latitude_fmt=":02.0f",
    map_extent=None,
    map_projection=ccrs.Mercator(),
    quiver_alpha=None,
    quiver_coarsen=10,
    quiver_color="k",
    quiver_visible=False,
    save_extra="",
    save_figure=False,
    save_path="",
    ship_color=None,
    ship_distance=None,
    ship_fade_concentric=True,
    ship_lon_lat=None,
    title="",
    vmin=None,
    vmax=None,
)
```

#### Loop through all days

To draw maps for every day of data in the dataset (each on a separate figure), use `surface_map_daily`:

```python
dc.plot.surface_map_daily(data, fvar, **kwargs)
```

The above automatically sets a consistent `vmin` and `vmax` across all the figures for each variable.  The `kwargs` are mostly the same as for `surface_map` above.

### Time series at a point

#### Show a single variable

`surface_timeseries` draws a time series with PCHIP interpolation between the points.  You first need to make the slice yourself so that `time` is the only remaining dimension in the dataset, as shown below.

```python
surface_point_slice = surface.sel(
    longitude=10.0, latitude=-38.4, method="nearest"
)  # for example

fig, ax = dc.plot.surface_timeseries(
    surface_point_slice,
    fvar,
    ax=None,
    dpi=150,
    figsize=[6.4, 4.8],
    draw_line=True,
    draw_points=True,
    interpolate_pchip=True,
    show_offset_text=True,
)
```

#### Grids with multiple variables

You can draw a grid of all variables of a particular type in a single figure:

```python
# Physical variables
fig, axs = dc.plot.surface_timeseries_grid_phys(
    surface_point_slice, dpi=150, figsize=[9.6, 7.2]
)

# Carbonate system variables
fig, axs = dc.plot.surface_timeseries_grid_co2(
    surface_point_slice, dpi=150, figsize=[9.6, 4.8]
)

# Biological variables
fig, axs = dc.plot.surface_timeseries_grid_bio(
    surface_point_slice, dpi=150, figsize=[9.6, 4.8]
)

# Nutrients
fig, axs = dc.plot.surface_timeseries_grid_nuts(
    surface_point_slice, dpi=150, figsize=[9.6, 4.8]
)
```

Or all the grids together (but there's currently no way to automatically save the figures that come out):

```python
dc.plot.surface_timeseries_grids(surface_point_slice, dpi=150)
```

#### Polar current plot

To draw a time series of current in polar coordinates:

```python
fig, ax = dc.plot.surface_currents(
    surface_point_slice, dpi=150, figsize=[6.4, 4.8]
)
```