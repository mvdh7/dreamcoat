from collections import namedtuple
import warnings
from matplotlib import pyplot as plt  # , patheffects as pe
import numpy as np
from scipy import interpolate
from sklearn import cluster
from . import meta, stats

Clustered = namedtuple("Clustered", ("x", "Y", "std", "std_unbiased", "count", "label"))


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
        # fontweight="bold",
        # path_effects=[pe.withStroke(linewidth=0.8, alpha=0.7, foreground="#ceb301")],
        ha="left",
        va="bottom",
        rotation=-90,
        transform=ax.transAxes,
        fontsize=6.5,
    )


def get_clusters(x, Y, cluster_bandwidth):
    """_summary_

    Parameters
    ----------
    x : np.ndarray
        _description_
    Y : np.ndarray
        _description_
    cluster_bandwitdh : float
        The ``bandwith`` kwarg for ``sklearn.cluster.MeanShift``.

    Returns
    -------
    np.ndarray
        _description_
    np.ndarray
        _description_
    np.ndarray
        _description_
    """
    # Eliminate NaNs from input data
    l = ~np.isnan(x)
    x = x[l]
    Y = Y[l]
    # Do MeanShift clustering of the x-variable
    ms = cluster.MeanShift(bandwidth=cluster_bandwidth)
    ms.fit(np.vstack(x))
    # Get y-variable values at x-variable cluster centres
    cn = ms.labels_
    cx = ms.cluster_centers_.ravel()
    try:
        cY = np.full((cx.size, Y.shape[1]), np.nan)
        cY_std = np.full((cx.size, Y.shape[1]), np.nan)
        cY_std_unbiased = np.full((cx.size, Y.shape[1]), np.nan)
        cY_count = np.full((cx.size, Y.shape[1]), np.nan)
    except IndexError:
        cY = np.full(cx.size, np.nan)
        cY_std = np.full(cx.size, np.nan)
        cY_std_unbiased = np.full(cx.size, np.nan)
        cY_count = np.full(cx.size, 0)
    for i in range(cx.size):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="Degrees of freedom <= 0 for slice."
            )
            warnings.filterwarnings("ignore", message="Mean of empty slice")
            cY[i] = np.nanmean(Y[cn == i], axis=0)
            cY_std[i] = np.nanstd(Y[cn == i], axis=0)
        cY_std_unbiased[i] = stats.std_unbiased(Y[cn == i], axis=0)
        try:
            cY_count[i] = np.sum(np.vstack(cn == i) & ~np.isnan(Y), axis=0)
        except ValueError:
            cY_count[i] = np.sum((cn == i) & ~np.isnan(Y), axis=0)
    # Sort clusters
    ci = np.argsort(cx)
    cx = cx[ci]
    cY = cY[ci]
    cY_std = cY_std[ci]
    cY_std_unbiased = cY_std_unbiased[ci]
    cY_count = cY_count[ci]
    return Clustered(cx, cY, cY_std, cY_std_unbiased, cY_count, cn)


def get_cluster_profile(x, y, cluster_bandwidth=5, linspace_num=100):
    """Get clusters from input data and interpolate them."""
    cl = get_clusters(x, y, cluster_bandwidth=cluster_bandwidth)
    # Generate PCHIP interpolation of clustered data
    px = np.linspace(np.min(cl.x), np.max(cl.x), num=linspace_num)
    py = interpolate.pchip(cl.x, cl.Y)(px)
    return px, py, cl.x, cl.Y, cl.label


def cluster_profile(
    x,
    y,
    data=None,
    ax=None,
    cluster_bandwidth=5,
    linspace_num=100,
    invert_xy=False,
    plot_kwargs=None,
    scatter_kwargs=None,
):
    """Get clusters from input data, interpolate them and plot the results.

    Parameters
    ----------
    x : np.ndarray
        The independent variable, which is used for clustering.
    y : np.ndarray
        The dependent variable.
    data : dict-like, optional
        If provided, ``x`` and ``y`` are treated as being the keys (so the independent
        variable data is in ``data[x]`` and dependent in ``data[y]``).
    ax : optional
        Matplotlib axes to draw the figure on.  If not provided, new axes are created.
    cluster_bandwitdh : float, optional
        The ``bandwith`` kwarg for ``sklearn.cluster.MeanShift``, by default 5.
    linspace_num : int, optional
        The ``num`` kwarg for ``np.linspace`` interpolation, by default 100.
    invert_xy : bool, optional
        Whether to plot the independent variable on the x-axis (False) or y-axis (True),
        by default False.
    plot_kwargs : dict, optional
        Any kwargs to pass to the ``ax.plot`` artist.
    scatter_kwargs : dict, optional
        Any kwargs to pass to the ``ax.scatter`` artist.
    """
    # Extract data from data, if required
    if data is not None:
        x = data[x]
        y = data[y]
    p_xy = get_cluster_profile(
        x,
        y,
        cluster_bandwidth=cluster_bandwidth,
        linspace_num=linspace_num,
    )[:2]
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()
    if plot_kwargs is None:
        plot_kwargs = {}
    if scatter_kwargs is None:
        scatter_kwargs = {}
    xy = (x, y)
    if invert_xy:
        p_xy = p_xy[::-1]
        xy = xy[::-1]
    ax.plot(*p_xy, **plot_kwargs)
    ax.scatter(*xy, **scatter_kwargs)
    add_credit(ax)
    return fig, ax
