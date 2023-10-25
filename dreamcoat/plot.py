import warnings
import numpy as np
from matplotlib import pyplot as plt, dates as mdates, patheffects as pe
from matplotlib.collections import LineCollection
from scipy.interpolate import pchip_interpolate
from cartopy import crs as ccrs, feature as cfeature
from cartopy.geodesic import Geodesic
from . import convert, meta


# TODO add EEZ and MPA boundaries to maps


def add_ship(
    ax,
    longitude,
    latitude,
    distance=None,
    color="xkcd:strawberry",
    fade_concentric=True,
    show_centre=True,
):
    """Add ship location and an optional circle at a certain distance around it.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
        The axes to add the credit to.
    longitude : float
        The ship's longitude in decimal degrees E.
    latitude : float
        The ship's latitude in decimal degrees N.
    distance : float, optional
        The distance around the ship to draw as a circle in km, by default None.  A list
        of multiple distances can be provided for a set of concentric circles.
    color : str, optional
        The colour in which to plot the ship and circle, by default "xkcd:strawberry".
    fade_concentric : bool, optional
        Whether to gradually fade out concentric circles if multiple distances are
        provided, by default True.
    show_centre : bool, optional
        Whether to show the central ship location, by default True.
    """
    if show_centre:
        ax.scatter(
            longitude,
            latitude,
            c=color,
            marker="+",
            lw=1,
            s=10,
            transform=ccrs.PlateCarree(),
        )
    if distance is not None:
        if np.isscalar(distance):
            distance = [distance]
        d_alpha = 1
        for d in distance:
            circle = Geodesic().circle(
                longitude, latitude, d * 1000, n_samples=180, endpoint=True
            )
            ax.plot(
                *circle.T, transform=ccrs.PlateCarree(), c=color, lw=1, alpha=d_alpha
            )
            if fade_concentric and np.size(distance) > 1:
                if np.size(distance) == 2:
                    d_alpha = 0.6
                else:
                    d_alpha -= 0.6 / (np.size(distance) - 1)


def add_credit(ax):
    """Add dreamcoat credit to axes.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
        The axes to add the credit to.
    """
    ax.text(
        1.008,
        0.01,
        "dreamcoat {}".format(meta.version_colour),
        # "dreamcoat {} (v{})".format(meta.version_colour, meta.version_number),
        c="xkcd:{}".format(meta.version_colour),
        fontweight="bold",
        path_effects=[pe.withStroke(linewidth=0.8, alpha=0.7, foreground="#ceb301")],
        ha="left",
        va="bottom",
        rotation=-90,
        transform=ax.transAxes,
        fontsize=6.5,
    )


styles = {
    "current_east": dict(
        cmap="RdBu",
        label="Eastwards current velocity / m/s",
        contrast="xkcd:kelly green",
    ),
    "current_north": dict(
        cmap="RdBu",
        label="Northwards current velocity / m/s",
        contrast="xkcd:kelly green",
    ),
    "current_speed": dict(
        cmap="cividis",
        label="Current speed / m/s",
        contrast="xkcd:strawberry",
    ),
    "current": dict(
        label="Current speed / m/s",
    ),
    "mld": dict(
        cmap="magma_r",
        label="Mixed layer depth / m",
        contrast="xkcd:kelly green",
    ),
    "salinity": dict(
        cmap="viridis",
        label="Practical salinity",
        contrast="xkcd:strawberry",
    ),
    "ssh": dict(
        cmap="BrBG_r",
        label="Sea surface height / m",
        contrast="xkcd:light purple",
    ),
    "theta": dict(
        cmap="plasma",
        label="Potential temperature / °C",
        contrast="xkcd:aqua",
    ),
    "talk": dict(
        cmap="viridis",
        label="Total alkalinity / µmol kg$^{-1}$",
        contrast="xkcd:strawberry",
    ),
    "dissic": dict(
        cmap="viridis",
        label="DIC / µmol kg$^{-1}$",
        contrast="xkcd:strawberry",
    ),
    "ph": dict(
        cmap="plasma_r",
        label="pH (total scale)",
        contrast="xkcd:cyan",
    ),
    "no3": dict(
        cmap="magma",
        label="Nitrate / µmol kg$^{-1}$",
        contrast="xkcd:cyan",
    ),
    "po4": dict(
        cmap="magma",
        label="Phosphate / µmol kg$^{-1}$",
        contrast="xkcd:cyan",
    ),
    "si": dict(
        cmap="magma",
        label="Silicate / µmol kg$^{-1}$",
        contrast="xkcd:cyan",
    ),
    "fe": dict(
        cmap="magma",
        label="Iron / nmol kg$^{-1}$",
        contrast="xkcd:cyan",
    ),
    "chl": dict(
        cmap="magma",
        label="Chlorophyll / mg m$^{-3}$",
        contrast="xkcd:cyan",
    ),
    "o2": dict(
        cmap="cividis",
        label="Oxygen / µmol kg$^{-1}$",
        contrast="xkcd:strawberry",
    ),
    "spco2": dict(
        cmap="viridis",
        label="Seawater $p$CO$_2$ / µatm",
        contrast="xkcd:strawberry",
    ),
    "nppv": dict(
        cmap="viridis",
        label="NPP / mg-C m$^{-3}$ day$^{-1}$",
        contrast="xkcd:strawberry",
    ),
    "phyc": dict(
        cmap="viridis",
        label="Phytoplankton / mmol-C m$^{-3}$",
        contrast="xkcd:strawberry",
    ),
    "density": dict(
        cmap="viridis",
        label="Density / kg dm$^{-3}$",
        contrast="xkcd:strawberry",
    ),
    "density_anomaly": dict(
        cmap="viridis",
        label="Density anomaly / kg m$^{-3}$",
        contrast="xkcd:strawberry",
    ),
    "aou": dict(
        cmap="RdBu",
        label="AOU / µmol kg$^{-1}$",
        contrast="xkcd:green",
    ),
    "pic": dict(
        cmap="viridis",
        label="PIC / mmol m$^{-3}$",
        contrast="xkcd:strawberry",
    ),
}


def _get_data_extent(data, map_extent):
    if map_extent is not None:
        data_extent = data.sel(
            longitude=slice(map_extent[0], map_extent[1]),
            latitude=slice(map_extent[2], map_extent[3]),
        )
        if data_extent[list(data_extent.variables)[0]].size == 0:
            print("map_extent falls outside data_extent - ignoring!")
            data_extent = data.copy()
    else:
        data_extent = data.copy()
    return data_extent


def _get_vmin_vmax(data_extent, color_zoom_factor):
    fstyles = {}
    # Variables centred on zero with vmin = -vmax
    for v in ["current_east", "current_north", "ssh", "aou"]:
        if v in data_extent:
            fstyles.update(
                {
                    v: dict(
                        vmin=-np.max(np.abs(data_extent[v])) * color_zoom_factor,
                        vmax=np.max(np.abs(data_extent[v])) * color_zoom_factor,
                    )
                }
            )
    # Variables with lowest values close to, but not lower than, zero
    for v in [
        "mld",
        "no3",
        "po4",
        "si",
        "fe",
        "chl",
        "nppv",
        "phyc",
        "current_speed",
        "pic",
    ]:
        if v in data_extent:
            fstyles.update(
                {
                    v: dict(
                        vmin=0,
                        vmax=data_extent[v].max() * color_zoom_factor,
                    )
                }
            )
    # Other variables
    for v in [
        "talk",
        "dissic",
        "ph",
        "o2",
        "spco2",
        "density",
        "theta",
        "salinity",
        "density_anomaly",
    ]:
        if v in data_extent:
            var_range = data_extent[v].max() - data_extent[v].min()
            fstyles.update(
                {
                    v: dict(
                        vmin=data_extent[v].min()
                        + var_range * (1 - color_zoom_factor) / 2,
                        vmax=data_extent[v].max()
                        - var_range * (1 - color_zoom_factor) / 2,
                    )
                }
            )
    return fstyles


def add_lon_lat_labels(
    ax, longitude_fmt=":03.0f", latitude_fmt=":02.0f", map_extent=None, fontsize=10
):
    if map_extent is not None:
        ax.set_extent(map_extent, crs=ccrs.PlateCarree())
    map_extent = ax.get_extent(crs=ccrs.PlateCarree())
    nsew = convert.extent_to_nsew(map_extent)
    ax.text(
        -0.003,
        0,
        ("{" + latitude_fmt + "}°{}").format(np.abs(map_extent[2]), nsew[2]),
        rotation=90,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=fontsize,
    )
    ax.text(
        -0.003,
        1,
        ("{" + latitude_fmt + "}°{}").format(np.abs(map_extent[3]), nsew[3]),
        rotation=90,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=fontsize,
    )
    ax.text(
        0,
        1.003,
        ("{" + longitude_fmt + "}°{}").format(np.abs(map_extent[0]), nsew[0]),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=fontsize,
    )
    ax.text(
        1,
        1.003,
        ("{" + longitude_fmt + "}°{}").format(np.abs(map_extent[1]), nsew[1]),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=fontsize,
    )
    return map_extent


def surface_map(
    data,
    fvar,
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
    contrast=None,
    ship_distance=None,
    ship_fade_concentric=True,
    ship_lon_lat=None,
    ship_show_centre=True,
    title="",
    vmin=None,
    vmax=None,
):
    """Plot a map of a surface physics or biogeochemistry dataset from CMEMS.

    Parameters
    ----------
    data : xarray Dataset.
        The dataset, opened with cmems.open_surface(), at a single time point.
    fvar : str
        Which variable from the dataset to plot.
    ax : matplotlib axes, optional
        An existing set of axes to plot onto, by default None.
    color_zoom_factor : float, optional
        What fraction of the total range of values to show on the colour bar, by default
        0.9.
    dpi : int, optional
        Figure resolution in dots per inch, by default 150.
    figsize : list, optional
        Figure size in inches, by default [6.4, 4.8].
    land_visible : bool, optional
        Whether to draw land areas from Natural Earth Data, by default True.
    longitude_fmt : str, optional
        Format of the longitude axis labels, by default ':03.0f'.
    latitude_fmt : str, optional
        Format of the latitude axis labels, by default ':02.0f'.
    map_extent : list, optional
        Map extents in decimal degrees, if different from those automatically generated
        by the dataset, by default None.
    map_projection : cartopy.crs, optional
        Cartopy projection system to use, by default ccrs.Mercator().
    quiver_alpha : float, optional
        Opacity of quiver arrows, by default None.
    quiver_coarsen : int, optional
        Coarsening factor for quiver arrows, by default 10.
    quiver_color : str, optional
        Colour of quiver arrows, by default "k".
    quiver_visible : bool, optional
        Whether to show the quiver arrows (they are always shown when plotting the
        'current_speed'), by default False.
    save_extra : str, optional
        Any extra text to append to the figure file name, by default "".
    save_figure : bool, optional
        Whether to save the figure to file, by default False.
    save_path : str, optional
        The file path to save the figure to, by default "".
    contrast : str, optional
        What colour to plot the ship's location in, by default None.
    ship_distance : float, optional
        Distance of the concentric circle lines around the ship's location in km, by
        default None.
    ship_fade_concentric : bool, optional
        Whether to fade out the concentric circle lines, by default True.
    ship_lon_lat : float, optional
        Longitude and latitude of the ship in decimal degrees N or E, by default None.
    ship_show_centre : bool, optional
        Whether to show the central ship location, by default True.
    title : str, optional
        Text for the figure title, by default "".
    vmin : float, optional
        Minimum value for the colour bar, by default None.
    vmax : float, optional
        Maximum value for the colour bar, by default None.

    Returns
    -------
    matplotlib figure
        The generated matplotlib figure.
    matplotlib axis
        The generated matplotlib axis.
    """
    # Set up dict of settings for all variables
    data_extent = _get_data_extent(data, map_extent)
    fstyles = _get_vmin_vmax(data_extent, color_zoom_factor)
    for k in fstyles:
        fstyles[k].update(styles[k])
    # Finalise settings for the selected variable
    fs = fstyles[fvar].copy()
    if contrast is not None:
        fs["contrast"] = contrast
    if vmin is not None:
        fs["vmin"] = vmin
    if vmax is not None:
        fs["vmax"] = vmax

    # Initialise the figure, if necessary
    if ax is None:
        fig = plt.figure(dpi=dpi, figsize=figsize)
        ax = fig.add_subplot(projection=map_projection)
    else:
        fig = ax.get_figure()

    # Plot the selected variable
    dplot = data[fvar].plot(
        ax=ax,
        add_colorbar=False,
        cmap=fs["cmap"],
        transform=ccrs.PlateCarree(),
        vmin=fs["vmin"],
        vmax=fs["vmax"],
    )

    # Determine how to extend the colorbar to show values outside its range
    extend = "neither"
    extend_sum = 0
    if data[fvar].min() < fs["vmin"]:
        extend = "min"
        extend_sum += 1
    if data[fvar].max() > fs["vmax"]:
        extend = "max"
        extend_sum += 1
    if extend_sum == 2:
        extend = "both"
    # Draw the colorbar
    cb = fig.colorbar(
        dplot,
        aspect=25,
        extend=extend,
        label=fstyles[fvar]["label"],
        location="bottom",
        pad=0.05,
        shrink=0.8,
    )
    # Draw quivers if showing current speed
    if fvar == "current_speed" or quiver_visible:
        # Unpack coarsening factor for quivers
        if np.size(quiver_coarsen) == 2:
            coarsen_x, coarsen_y = quiver_coarsen
        else:
            coarsen_x = coarsen_y = quiver_coarsen
        # Coarsen the quiver data
        data_coarse = data_extent.coarsen(
            longitude=coarsen_x, latitude=coarsen_y, boundary="trim"
        ).mean()
        # Determine opacity of the quivers, if not provided
        if quiver_alpha is None:
            if quiver_color == "w":
                quiver_alpha = 0.25
            elif quiver_color == "k":
                quiver_alpha = 0.6
            else:
                quiver_alpha = 0.5
        # Draw quiver plot
        lon, lat = np.meshgrid(data_coarse.longitude, data_coarse.latitude)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            ax.quiver(
                lon,
                lat,
                data_coarse.current_east.data,
                data_coarse.current_north.data,
                alpha=quiver_alpha,
                color=quiver_color,
                transform=ccrs.PlateCarree(),
            )

    # Show land areas
    if land_visible:
        ax.add_feature(
            cfeature.NaturalEarthFeature("physical", "land", "10m"),
            edgecolor="none",
            facecolor=np.array([1, 1, 1]) * 0.1,
        )

    # Add latitude / longitude extents and labels
    map_extent = add_lon_lat_labels(
        ax,
        longitude_fmt=longitude_fmt,
        latitude_fmt=latitude_fmt,
        map_extent=map_extent,
    )

    # Misc. additions and settings
    ax.set_title(title)
    if ship_lon_lat is not None:
        add_ship(
            ax,
            *ship_lon_lat,
            distance=ship_distance,
            color=fs["contrast"],
            fade_concentric=ship_fade_concentric,
            show_centre=ship_show_centre,
        )
    add_credit(ax)
    ax.set_extent(map_extent, crs=ccrs.PlateCarree())  # in case it's changed
    fig.tight_layout()

    # Save to file, if requested, and finish
    if save_figure:
        fig.savefig(save_path + "surface_{}{}.png".format(fvar, save_extra))
    return fig, ax


def surface_map_daily(
    data,
    fvar,
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
    save_figure=False,
    save_path="",
    contrast=None,
    ship_distance=None,
    ship_fade_concentric=True,
    ship_lon_lat=None,
    vmin=None,
    vmax=None,
):
    """Plot daily maps of a surface physics or biogeochemistry dataset from CMEMS.

    Parameters
    ----------
    data : xarray Dataset.
        The dataset, opened with cmems.open_surface().
    fvar : str
        Which variable from the dataset to plot.
    color_zoom_factor : float, optional
        What fraction of the total range of values to show on the colour bar, by default
        0.9.
    dpi : int, optional
        Figure resolution in dots per inch, by default 150.
    figsize : list, optional
        Figure size in inches, by default [6.4, 4.8].
    land_visible : bool, optional
        Whether to draw land areas from Natural Earth Data, by default True.
    longitude_fmt : str, optional
        Format of the longitude axis labels, by default ':03.0f'.
    latitude_fmt : str, optional
        Format of the latitude axis labels, by default ':02.0f'.
    map_extent : list, optional
        Map extents in decimal degrees, if different from those automatically generated
        by the dataset, by default None.
    map_projection : cartopy.crs, optional
        Cartopy projection system to use, by default ccrs.Mercator().
    quiver_alpha : float, optional
        Opacity of quiver arrows, by default None.
    quiver_coarsen : int, optional
        Coarsening factor for quiver arrows, by default 10.
    quiver_color : str, optional
        Colour of quiver arrows, by default "k".
    quiver_visible : bool, optional
        Whether to show the quiver arrows (they are always shown when plotting the
        'current_speed'), by default False.
    save_figure : bool, optional
        Whether to save the figure to file, by default False.
    save_path : str, optional
        The file path to save the figure to, by default "".
    contrast : str, optional
        What colour to plot the ship's location in, by default None.
    ship_distance : float, optional
        Distance of the concentric circle lines around the ship's location in km, by
        default None.
    ship_fade_concentric : bool, optional
        Whether to fade out the concentric circle lines, by default True.
    ship_lon_lat : float, optional
        Longitude and latitude of the ship in decimal degrees N or E, by default None.
    vmin : float, optional
        Minimum value for the colour bar, by default None.
    vmax : float, optional
        Maximum value for the colour bar, by default None.
    """
    # Get consistent vmin and vmax to use across all figures
    data_extent = _get_data_extent(data, map_extent)
    fstyles = _get_vmin_vmax(data_extent, color_zoom_factor)
    if vmin is None:
        vmin = fstyles[fvar]["vmin"]
    if vmax is None:
        vmax = fstyles[fvar]["vmax"]
    # Loop through days and make plots
    for i in range(data.time.size):
        _, ax = surface_map(
            data_extent.isel(time=i),
            fvar,
            color_zoom_factor=color_zoom_factor,
            dpi=dpi,
            figsize=figsize,
            land_visible=land_visible,
            longitude_fmt=longitude_fmt,
            latitude_fmt=latitude_fmt,
            map_extent=None,  # already handled above
            map_projection=map_projection,
            quiver_alpha=quiver_alpha,
            quiver_coarsen=quiver_coarsen,
            quiver_color=quiver_color,
            quiver_visible=quiver_visible,
            save_extra="_{:02.0f}".format(i),
            save_figure=save_figure,
            save_path=save_path,
            contrast=contrast,
            ship_distance=ship_distance,
            ship_fade_concentric=ship_fade_concentric,
            ship_lon_lat=ship_lon_lat,
            title=str(data.time.data[i]).split("T")[0],
            vmin=vmin,
            vmax=vmax,
        )
        plt.show()
        plt.close()


def _get_surface_timeseries_line(data, fvar, interpolate_pchip):
    if interpolate_pchip:
        ndates = mdates.date2num(data.time)
        ix = np.linspace(ndates[0], ndates[-1], num=ndates.size * 10)
        fx = mdates.num2date(ix)
        fy = pchip_interpolate(ndates, data[fvar].data, ix)
    else:
        fx = data.time.data
        fy = data[fvar].data
    return fx, fy


def surface_timeseries(
    data,
    fvar,
    ax=None,
    dpi=150,
    figsize=[6.4, 4.8],
    draw_line=True,
    draw_points=True,
    interpolate_pchip=True,
    show_offset_text=True,
):
    """Show a timeseries of data from a given location.

    Parameters
    ----------
    data : xarray Dataset.
        The dataset, opened with cmems.open_surface(), at a single time point.
    fvar : str
        Which variable from the dataset to plot.
    ax : matplotlib axes, optional
        An existing set of axes to plot onto, by default None.
    dpi : int, optional
        Figure resolution in dots per inch, by default 150.
    figsize : list, optional
        Figure size in inches, by default [6.4, 4.8].
    draw_line : bool, optional
        Whether to draw lines between the data points, by default True.
    draw_points : bool, optional
        Whether to show the data points with markers, by default True.
    interpolate_pchip : bool, optional
        Whether to use PCHIP interpolation for the lines, by default True.
    show_offset_text : bool, optional
        Whether to show the offset text on the date axis, by default True.

    Returns
    -------
    matplotlib figure
        The generated matplotlib figure.
    matplotlib axis
        The generated matplotlib axis.
    """
    # Create figure, if necessary
    if ax is None:
        fig, ax = plt.subplots(dpi=dpi, figsize=figsize)
    else:
        fig = ax.get_figure()

    # Plot data
    if draw_line:
        if fvar == "current":
            fx, fy = _get_surface_timeseries_line(
                data, "current_north", interpolate_pchip
            )
            ax.plot(fx, fy, c="xkcd:cerulean", alpha=0.8, label="North")
            fx, fy = _get_surface_timeseries_line(
                data, "current_east", interpolate_pchip
            )
            ax.plot(fx, fy, c="xkcd:tangerine", alpha=0.8, label="East")
            fx, fy = _get_surface_timeseries_line(
                data, "current_speed", interpolate_pchip
            )
            ax.plot(fx, fy, c="xkcd:almost black", alpha=0.8, label="Total", zorder=10)
            ax.legend()
        else:
            fx, fy = _get_surface_timeseries_line(data, fvar, interpolate_pchip)
            ax.plot(fx, fy, c="xkcd:navy")
    if draw_points:
        if fvar == "current":
            ax.scatter(
                "time",
                "current_north",
                alpha=0.8,
                c="xkcd:cerulean",
                data=data,
                label=None,
                s=20,
            )
            ax.scatter(
                "time",
                "current_east",
                alpha=0.8,
                c="xkcd:tangerine",
                data=data,
                label=None,
                s=20,
            )
            ax.scatter(
                "time",
                "current_speed",
                alpha=0.8,
                c="xkcd:almost black",
                data=data,
                label=None,
                s=20,
                zorder=10,
            )
        else:
            ax.scatter("time", fvar, data=data, c="xkcd:navy", s=20)
    if fvar == "mld":
        ax.invert_yaxis()

    # Add x-axis zero line if needed
    if fvar == "current" or (
        fvar in ["current_north", "current_east", "ssh"]
        and data[fvar].max() > 0
        and data[fvar].min() < 0
    ):
        ax.axhline(0, c="k", lw=0.8)

    # Axis display settings
    ax.set_ylabel(styles[fvar]["label"])
    locator = mdates.AutoDateLocator(minticks=5, maxticks=12)
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    if not show_offset_text:
        ax.xaxis.get_offset_text().set_visible(False)
    add_credit(ax)

    return fig, ax


def surface_timeseries_grid_phys(data, dpi=150, figsize=[9.6, 7.2]):
    """Draw time-series plots for all physical variables as subplots of a single figure.

    Parameters
    ----------
    data : xarray Dataset.
        The dataset, opened with cmems.open_surface(), at a single point in space.
    dpi : int, optional
        Figure resolution in dots per inch, by default 150.
    figsize : list, optional
        Figure size in inches, by default [9.6, 7.2].

    Returns
    -------
    matplotlib figure
        The generated matplotlib figure.
    matplotlib axis
        The generated matplotlib axes.
    """
    fig, axs = plt.subplots(nrows=3, ncols=2, dpi=dpi, figsize=figsize)
    fvars = ["theta", "salinity", "ssh", "density_anomaly", "mld", "current"]
    letters = "abcdef"
    for ax, fvar, letter in zip(axs.ravel(), fvars, letters):
        if fvar:
            surface_timeseries(data, fvar, ax=ax)
            ax.text(0, 1.05, "(" + letter + ")", transform=ax.transAxes)
        else:
            ax.set_visible(False)
    fig.tight_layout()
    return fig, axs


def surface_timeseries_grid_co2(data, dpi=150, figsize=[9.6, 4.8]):
    """Draw time-series plots for all marine carbonate system variables as subplots of a
    single figure.

    Parameters
    ----------
    data : xarray Dataset.
        The dataset, opened with cmems.open_surface(), at a single point in space.
    dpi : int, optional
        Figure resolution in dots per inch, by default 150.
    figsize : list, optional
        Figure size in inches, by default [9.6, 7.2].

    Returns
    -------
    matplotlib figure
        The generated matplotlib figure.
    matplotlib axis
        The generated matplotlib axes.
    """
    fig, axs = plt.subplots(nrows=2, ncols=2, dpi=dpi, figsize=figsize)
    fvars = ["talk", "dissic", "ph", "spco2"]
    letters = "abcdef"
    for ax, fvar, letter in zip(axs.ravel(), fvars, letters):
        if fvar:
            surface_timeseries(data, fvar, ax=ax)
            ax.text(0, 1.05, "(" + letter + ")", transform=ax.transAxes)
        else:
            ax.set_visible(False)
    fig.tight_layout()
    return fig, axs


def surface_timeseries_grid_bio(data, dpi=150, figsize=[9.6, 4.8]):
    """Draw time-series plots for all biological variables as subplots of a single
    figure.

    Parameters
    ----------
    data : xarray Dataset.
        The dataset, opened with cmems.open_surface(), at a single point in space.
    dpi : int, optional
        Figure resolution in dots per inch, by default 150.
    figsize : list, optional
        Figure size in inches, by default [9.6, 7.2].

    Returns
    -------
    matplotlib figure
        The generated matplotlib figure.
    matplotlib axis
        The generated matplotlib axes.
    """
    fig, axs = plt.subplots(nrows=2, ncols=2, dpi=dpi, figsize=figsize)
    fvars = ["phyc", "nppv", "chl", "o2"]
    letters = "abcdef"
    for ax, fvar, letter in zip(axs.ravel(), fvars, letters):
        if fvar:
            surface_timeseries(data, fvar, ax=ax)
            ax.text(0, 1.05, "(" + letter + ")", transform=ax.transAxes)
        else:
            ax.set_visible(False)
    fig.tight_layout()
    return fig, axs


def surface_timeseries_grid_nuts(data, dpi=150, figsize=[9.6, 4.8]):
    """Draw time-series plots for all nutrient variables as subplots of a single figure.

    Parameters
    ----------
    data : xarray Dataset.
        The dataset, opened with cmems.open_surface(), at a single point in space.
    dpi : int, optional
        Figure resolution in dots per inch, by default 150.
    figsize : list, optional
        Figure size in inches, by default [9.6, 7.2].

    Returns
    -------
    matplotlib figure
        The generated matplotlib figure.
    matplotlib axis
        The generated matplotlib axes.
    """
    fig, axs = plt.subplots(nrows=2, ncols=2, dpi=dpi, figsize=figsize)
    fvars = ["no3", "po4", "si", "fe"]
    letters = "abcdef"
    for ax, fvar, letter in zip(axs.ravel(), fvars, letters):
        if fvar:
            surface_timeseries(data, fvar, ax=ax)
            ax.text(0, 1.05, "(" + letter + ")", transform=ax.transAxes)
        else:
            ax.set_visible(False)
    fig.tight_layout()
    return fig, axs


def surface_timeseries_grids(data, dpi=150):
    """Draw all time-series grid figures.

    Parameters
    ----------
    data : xarray Dataset.
        The dataset, opened with cmems.open_surface(), at a single point in space.
    dpi : int, optional
        Figure resolution in dots per inch, by default 150.
    """
    surface_timeseries_grid_phys(data, dpi=dpi)
    surface_timeseries_grid_co2(data, dpi=dpi)
    surface_timeseries_grid_bio(data, dpi=dpi)
    surface_timeseries_grid_nuts(data, dpi=dpi)


def _get_surface_currents_line(data, interpolate_pchip):
    if interpolate_pchip:
        ndates = mdates.date2num(data.time)
        it = np.linspace(ndates[0], ndates[-1], num=ndates.size * 10)
        # ft = mdates.num2date(it)
        fx = pchip_interpolate(ndates, data.current_east.data, it)
        fy = pchip_interpolate(ndates, data.current_north.data, it)
    else:
        it = mdates.date2num(data.time)
        fx = data.current_east.data
        fy = data.current_north.data
    ftheta, frho = convert.cartesian_to_polar(fx, fy)
    return it, ftheta, frho


def surface_currents(data, dpi=150, figsize=[6.4, 4.8]):
    """Draw a polar plot of surface currents, coloured by date.

    Parameters
    ----------
    data : xarray Dataset.
        The dataset, opened with cmems.open_surface(), at a single time point.
    dpi : int, optional
        Figure resolution in dots per inch, by default 150.
    figsize : list, optional
        Figure size in inches, by default [6.4, 4.8].

    Returns
    -------
    matplotlib figure
        The generated matplotlib figure.
    matplotlib axis
        The generated matplotlib axis.
    """
    # Initialise figure
    fig = plt.figure(dpi=150)
    ax = fig.add_subplot(111, projection="polar")

    # Get points in cartesian coordinates
    it, ftheta, frho = _get_surface_currents_line(data, True)
    points = np.array([ftheta, frho]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Correct angle orientation for plot (north up and clockwise increase)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    # Create and draw coloured lines
    norm = plt.Normalize(it.min(), it.max())
    lc = LineCollection(
        segments,
        cmap="viridis",
        norm=norm,
        zorder=11,
        path_effects=[pe.Stroke(capstyle="round")],
    )
    lc.set_array(it)
    lc.set_linewidth(3)
    line = ax.add_collection(lc)

    # Scatter points
    theta, rho = convert.cartesian_to_polar(
        data.current_east.data, data.current_north.data
    )
    ax.scatter(theta, rho, c=mdates.date2num(data.time.data), s=25, zorder=10)

    # # Colorbar
    # plt.colorbar(line, location='right')

    # Finish off
    add_credit(ax)
    ax.set_ylim([0, np.max(frho * 1.1)])
    fig.tight_layout()
    return fig, ax
