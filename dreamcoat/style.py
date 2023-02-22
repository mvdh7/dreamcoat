import numpy as np


def default_if_None(value, default):
    if value is None:
        return default
    else:
        return value


contrasts = {
    "viridis": "xkcd:strawberry",
    "RdBu": "xkcd:kelly green",
    "cividis": "xkcd:strawberry",
    "magma": "xkcd:kelly green",
    "BrBG": "xkcd:light purple",
    "plasma": "xkcd:aqua",
}
contrasts.update({k + "_r": v for k, v in contrasts.items()})


def get_contrast(cmap):
    """Return an appropriate contrast colour for each cmap."""
    if cmap in contrasts:
        return contrasts[cmap]
    else:
        return "xkcd:strawberry"


class Style:
    def __init__(
        self,
        alpha=0.9,
        cmap="viridis",
        contrast=None,
        label=None,
        name=None,
    ):
        self.alpha = alpha
        self.cmap = cmap
        if contrast is None:
            contrast = get_contrast(cmap)
        self.contrast = contrast
        self.label = label
        self.name = name

    def __repr__(self):
        return """dreamcoat.Style instance
  alpha    = {alpha}
  cmap     = {cmap}
  contrast = {contrast}
  label    = {label}""".format(
            alpha=self.alpha, cmap=self.cmap, contrast=self.contrast, label=self.label
        )

    def get_vmin_vmax(self, values, percentile=1):
        """Compute appropriate vmin and vmax values of a given variable for given values."""
        values = np.array(values)
        values = values[~np.isnan(values)]
        # Variables centred on zero with vmin = -vmax
        if self.name in ["current_east", "current_north", "ssh", "aou"]:
            lower = np.percentile(values, percentile)
            upper = np.percentile(values, 100 - percentile)
            limit = np.max((np.abs(lower), np.abs(upper)))
            vmin = -limit
            vmax = limit
        # Variables with lowest values close to, but not lower than, zero
        elif self.name in [
            "mld",
            "nitrate",
            "phosphate",
            "silicate",
            "iron",
            "chlorophyll",
            "npp",
            "phytoplankton",
            "current_speed",
            "pic",
        ]:
            vmin = 0
            vmax = np.percentile(values, 100 - percentile)
        # Other variables
        elif self.name in [
            "alkalinity",
            "dic",
            "pH",
            "oxygen",
            "pco2",
            "density",
            "theta",
            "salinity",
            "density_anomaly",
        ]:
            vmin = np.percentile(values, percentile)
            vmax = np.percentile(values, 100 - percentile)
        # np.percentile(modis.pic.data[~np.isnan(modis.pic.data)], 99.5)
        return dict(vmin=vmin, vmax=vmax)


# Create a Style for each variable
default = Style()
current_east = Style(
    name="current_east",
    cmap="RdBu",
    label="Eastwards current velocity / m/s",
)
current_north = Style(
    name="current_north",
    cmap="RdBu",
    label="Northwards current velocity / m/s",
)
current_speed = Style(
    name="current_speed",
    cmap="cividis",
    label="Current speed / m/s",
)
current = Style(
    name="current",
    label="Current speed / m/s",
)
depth = Style(name="depth", label="Depth / m")
mld = Style(
    name="mld",
    cmap="magma_r",
    label="Mixed layer depth / m",
)
salinity = Style(
    name="salinity",
    label="Practical salinity",
)
ssh = Style(
    name="ssh",
    cmap="BrBG_r",
    label="Sea surface height / m",
)
theta = Style(
    name="theta",
    cmap="plasma",
    label="Potential temperature / °C",
)
temperature = Style(
    name="temperature",
    cmap="plasma",
    label="Temperature / °C",
)
alkalinity = Style(
    name="alkalinity",
    label="Total alkalinity / µmol kg$^{-1}$",
)
dic = Style(
    name="dic",
    label="DIC / µmol kg$^{-1}$",
)
pH = Style(
    name="pH",
    cmap="plasma_r",
    label="pH (total scale)",
)
nitrate = Style(
    name="nitrate",
    cmap="magma",
    label="Nitrate / µmol kg$^{-1}$",
)
phosphate = Style(
    name="phosphate",
    cmap="magma",
    label="Phosphate / µmol kg$^{-1}$",
)
silicate = Style(
    name="silicate",
    cmap="magma",
    label="Silicate / µmol kg$^{-1}$",
)
iron = Style(
    name="iron",
    cmap="magma",
    label="Iron / nmol kg$^{-1}$",
)
chlorophyll = Style(
    name="chlorophyll",
    label="Chlorophyll / mg m$^{-3}$",
)
oxygen = Style(
    name="oxygen",
    cmap="cividis",
    label="Oxygen / µmol kg$^{-1}$",
)
pCO2 = Style(
    name="pCO2",
    label="Seawater $p$CO$_2$ / µatm",
)
npp = Style(
    name="npp",
    label="NPP / mg-C m$^{-3}$ day$^{-1}$",
)
phytoplankton = Style(
    name="phytoplankton",
    label="Phytoplankton / mmol-C m$^{-3}$",
)
density = Style(
    name="density",
    label="Density / kg dm$^{-3}$",
)
density_anomaly = Style(
    name="density_anomaly",
    label="Density anomaly / kg m$^{-3}$",
)
aou = Style(
    name="aou",
    cmap="RdBu",
    label="AOU / µmol kg$^{-1}$",
)
pic = Style(
    name="pic",
    label="PIC / mmol m$^{-3}$",
)

# Assemble the Styles into a dict for convenience
styles = {
    "default": default,
    "depth": depth,
    "dic": dic,
    "theta": theta,
    "current_east": current_east,
    "current_north": current_north,
    "current_speed": current_speed,
    "current": current,
    "mld": mld,
    "salinity": salinity,
    "ssh": ssh,
    "theta": theta,
    "alkalinity": alkalinity,
    "pH": pH,
    "nitrate": nitrate,
    "phosphate": phosphate,
    "silicate": silicate,
    "iron": iron,
    "chlorophyll": chlorophyll,
    "oxygen": oxygen,
    "pCO2": pCO2,
    "npp": npp,
    "phytoplankton": phytoplankton,
    "density": density,
    "density_anomaly": density_anomaly,
    "aou": aou,
    "pic": pic,
    "temperature": temperature,
}
