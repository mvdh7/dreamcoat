import warnings
import numpy as np
from matplotlib import pyplot as plt
from cartopy import crs as ccrs, feature as cfeature
from cartopy.geodesic import Geodesic
from . import convert, meta


def add_ship(
    ax,
    longitude,
    latitude,
    distance=None,
    color="xkcd:strawberry",
    fade_concentric=True,
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
    """
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
        1.005,
        0,
        "dreamcoat v{}".format(meta.version_number),
        c="xkcd:{}".format(meta.version_colour),
        ha="left",
        va="bottom",
        rotation=-90,
        transform=ax.transAxes,
        fontsize=7,
    )


def surphys(
    data,
    fvar,
    itime=0,
    color_zoom_factor=0.9,
    dpi=300,
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
    ship_color=None,
    ship_distance=None,
    ship_fade_concentric=True,
    ship_latitude=None,
    ship_longitude=None,
    vmin=None,
    vmax=None,
):
    """Plot a surface physics dataset from CMEMS.

    Parameters
    ----------
    data : xarray Dataset.
        The dataset, opened with cmems.open_surphys().
    fvar : str
        Which variable from the dataset to plot: 'theta', 'salinity', 'mld', 'ssh',
        'current_east', 'current_north' or 'current_speed'.
    itime : int, optional
        Which timestep from the dataset to plot, by default 0.
    color_zoom_factor : float, optional
        What fraction of the total range of values to show on the colour bar, by default
        0.9.
    dpi : int, optional
        Figure resolution in dots per inch, by default 300.
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
    ship_color : str, optional
        What colour to plot the ship's location in, by default None.
    ship_distance : float, optional
        Distance of the concentric circle lines around the ship's location in km, by default None.
    ship_fade_concentric : bool, optional
        Whether to fade out the concentric circle lines, by default True.
    ship_latitude : float, optional
        Latitude of the ship in decimal degrees N, by default None.
    ship_longitude : float, optional
        Longitude of the ship in decimal degrees E, by default None.
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
    salinity_range = data.salinity.max() - data.salinity.min()
    theta_range = data.theta.max() - data.salinity.min()
    fvar_settings = {
        "current_east": dict(
            cmap="RdBu",
            label="Eastwards current velocity / m/s",
            ship_color="xkcd:kelly green",
            vmin=-np.max(np.abs(data.current_east)) * color_zoom_factor,
            vmax=np.max(np.abs(data.current_east)) * color_zoom_factor,
        ),
        "current_north": dict(
            cmap="RdBu",
            label="Northwards current velocity / m/s",
            ship_color="xkcd:kelly green",
            vmin=-np.max(np.abs(data.current_north)) * color_zoom_factor,
            vmax=np.max(np.abs(data.current_north)) * color_zoom_factor,
        ),
        "current_speed": dict(
            cmap="cividis",
            label="Current speed / m/s",
            ship_color="xkcd:strawberry",
            vmin=0,
            vmax=data.current_speed.max() * color_zoom_factor,
        ),
        "mld": dict(
            cmap="magma_r",
            label="Mixed layer depth / m",
            ship_color="xkcd:kelly green",
            vmin=0,
            vmax=data.mld.max() * color_zoom_factor,
        ),
        "salinity": dict(
            cmap="viridis",
            label="Practical salinity",
            ship_color="xkcd:strawberry",
            vmin=data.salinity.min() + salinity_range * (1 - color_zoom_factor) / 2,
            vmax=data.salinity.max() - salinity_range * (1 - color_zoom_factor) / 2,
        ),
        "ssh": dict(
            cmap="BrBG_r",
            label="Sea surface height / m",
            ship_color="xkcd:light purple",
            vmin=-np.max(np.abs(data.ssh)) * color_zoom_factor,
            vmax=np.max(np.abs(data.ssh)) * color_zoom_factor,
        ),
        "theta": dict(
            cmap="plasma",
            label="Potential temperature / °C",
            ship_color="xkcd:aqua",
            vmin=data.theta.min() + theta_range * (1 - color_zoom_factor) / 2,
            vmax=data.theta.max() - theta_range * (1 - color_zoom_factor) / 2,
        ),
    }

    # Finalise settings for the selected variable
    fs = fvar_settings[fvar].copy()
    if ship_color is not None:
        fs["ship_color"] = ship_color
    if vmin is not None:
        fs["vmin"] = vmin
    if vmax is not None:
        fs["vmax"] = vmax

    # Initialise the figure
    fig = plt.figure(dpi=dpi, figsize=figsize)
    ax = fig.add_subplot(projection=map_projection)

    # Plot the selected variable
    dplot = (
        data[fvar]
        .isel(time=itime)
        .plot(
            ax=ax,
            add_colorbar=False,
            cmap=fs["cmap"],
            transform=ccrs.PlateCarree(),
            vmin=fs["vmin"],
            vmax=fs["vmax"],
        )
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
    plt.colorbar(
        dplot,
        aspect=25,
        extend=extend,
        label=fvar_settings[fvar]["label"],
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
        data_coarse = (
            data.isel(time=itime, depth=0)
            .coarsen(longitude=coarsen_x, latitude=coarsen_y, boundary="trim")
            .mean()
        )
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
            cfeature.NaturalEarthFeature("dataical", "land", "10m"),
            edgecolor="none",
            facecolor=np.array([1, 1, 1]) * 0.1,
        )

    # Add latitude / longitude extents and labels
    if map_extent is not None:
        ax.set_extent(map_extent, crs=ccrs.PlateCarree())
    map_extent = ax.get_extent(crs=ccrs.PlateCarree())
    nsew = convert.extent_to_nsew(map_extent)
    ax.text(
        -0.003,
        0,
        ("{" + longitude_fmt + "}°{}").format(np.abs(map_extent[2]), nsew[2]),
        rotation=90,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
    )
    ax.text(
        -0.003,
        1,
        ("{" + longitude_fmt + "}°{}").format(np.abs(map_extent[3]), nsew[3]),
        rotation=90,
        transform=ax.transAxes,
        ha="right",
        va="top",
    )
    ax.text(
        0,
        1.003,
        ("{" + latitude_fmt + "}°{}").format(np.abs(map_extent[0]), nsew[0]),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
    )
    ax.text(
        1,
        1.003,
        ("{" + latitude_fmt + "}°{}").format(np.abs(map_extent[1]), nsew[1]),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
    )

    # Misc. additions and settings
    ax.set_title(np.datetime_as_string(data.time[itime].data, "D"))
    add_ship(
        ax,
        ship_longitude,
        ship_latitude,
        distance=ship_distance,
        color=fs["ship_color"],
        fade_concentric=ship_fade_concentric,
    )
    add_credit(ax)
    ax.set_extent(map_extent, crs=ccrs.PlateCarree())  # in case it's changed
    plt.tight_layout()

    # Save to file, if requested
    if save_figure:
        plt.savefig(
            save_path
            + "surdata_{}_{}.png".format(
                fvar, np.datetime_as_string(data.time[itime].data, "D")
            )
        )

    return fig, ax
