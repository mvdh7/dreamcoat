from matplotlib import pyplot as plt, patheffects as pe
from . import meta


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
