import pandas as pd
import numpy as np
import networkx as nx
from scipy import interpolate
from matplotlib import pyplot as plt
from cartopy import crs as ccrs, feature as cfeature
from . import maps, plot


class CruiseGraph(nx.Graph):
    """A graph representation of a cruise for computing neutral surfaces.

    Parameters
    ----------
    ctdz : pd.DataFrame
        CTD sensor data for a cruise binned at regular depth or pressure intervals and
        including only downcast data.  Must contain at least the following columns:
            station : int or str
                An identifier that is unique for each CTD station.
            latitude, longitude : float
                Latitude and longitude in decimal degrees (N and E positive).
            pressure : float
                Hydrostatic pressure in dbar.
            theta : float
                Potential temperature in °C.
            salinity : float
                Practical salinity.
        The ctdz dataframe is stored in CruiseGraph.ctdz.
    edges : array-like, optional
        List of adjacency edges between pairs of stations, in the format required by
        nx.Graph.add_edges_from(), where edges here = ebunch_to_add there.

    Methods
    -------
    get_surfaces
        Compute all omega surfaces relative to a given reference station.
    get_surfaces_all
        Compute all omega surfaces with every station used as the reference separately.
    plot_surface_map
        Draw the adjacency graph and pressure on a given surface.
    """

    def __init__(self, ctdz, edges=None):
        super().__init__()
        self.ctdz = ctdz.copy()
        self._get_veronis_here()
        self._get_stations()
        self.surfaces_raw = {}
        self.surfaces = {}
        if edges is not None:
            self.add_edges_from(edges)

    def _get_veronis_here(self):
        # Calculate Veronis densities with each station referenced to itself
        # from neutralocean.label import veronis  # for neutralocean v2.2.0
        from neutralocean.label import veronis_density  # for neutralocean v2.1.3

        self.ctdz["veronis_here"] = np.nan
        for s, station in self.ctdz.groupby("station"):
            for i, row in station.iterrows():
                # For neutralocean v2.1.3:
                self.ctdz.loc[i, "veronis_here"] = veronis_density(
                    station.salinity.values,
                    station.theta.values,
                    station.pressure.values,
                    row.pressure,
                    eos="jmdfwg06",
                )
                # # For neutralocean v2.2.0:
                # self.ctdz.loc[i, "veronis_here"] = veronis(
                #     row.pressure,
                #     station.salinity.values,
                #     station.theta.values,
                #     station.pressure.values,
                #     eos="jmdfwg06",
                # )
        self.ctdz["veronis_here"] -= 1000

    def _get_stations(self):
        self.stations = (
            self.ctdz[["station", "longitude", "latitude"]].groupby("station").mean()
        )
        self.stations["npts"] = (
            self.ctdz[["station", "longitude"]].groupby("station").count()
        )
        self.add_nodes_from(self.stations.index)

    def _get_edge_distances(self):
        distances = {}
        for _a, _b in self.edges:
            row_a = self.stations.loc[_a]
            row_b = self.stations.loc[_b]
            distances[(_a, _b)] = maps.get_distance_geodesic(
                (row_a.longitude, row_a.latitude),
                (row_b.longitude, row_b.latitude),
            )
        nx.set_edge_attributes(self, distances, "distance")

    def add_edges_from(self, ebunch_to_add, **attr):
        assert np.all(
            np.isin(ebunch_to_add, self.stations.index)
        ), "All nodes must be in the existing list of stations!"
        super().add_edges_from(ebunch_to_add, **attr)
        self._get_edge_distances()
        self._get_grid()
        self._get_stp()

    def add_edge(self, u_of_edge, v_of_edge, **attr):
        assert (
            u_of_edge in self.stations.index and v_of_edge in self.stations.index
        ), "Both u_of_edge and v_of_edge must be in the existing list of stations!"
        super().add_edge(u_of_edge, v_of_edge, **attr)
        self._get_edge_distances()
        self._get_grid()
        self._get_stp()

    def _get_grid(self):
        a = []
        b = []
        distances = []
        for edge in self.edges:
            _a, _b = edge
            a.append((self.stations.index == _a).nonzero()[0][0])
            b.append((self.stations.index == _b).nonzero()[0][0])
            distances.append(self.edges[edge]["distance"])
        self.grid = {
            "edges": np.array([a, b]),
            "dist": np.array(distances),
            "distperp": np.ones_like(distances),
        }

    def _get_stp(self):
        stp_vars = ["salinity", "theta", "pressure"]
        stp = {}
        for v in stp_vars:
            stp[v] = np.full((self.stations.shape[0], self.stations.npts.max()), np.nan)
            for i, s in enumerate(self.stations.index):
                S = self.ctdz.station == s
                stp[v][i, : S.sum()] = self.ctdz[v][S].values
        self.stp = [stp[v] for v in stp_vars]

    def get_surfaces(self, station_ref, cutoff=5):
        """Compute all omega surfaces relative to a given reference station.

        After computing, surfaces are processed to remove outliers and smoothed so that
        each surface increases monotonically.  Veronis density labels are assigned from
        the reference station.

        Raw computed surfaces are stored in CruiseGraph.surfaces_raw[station_ref], and
        processed surfacs are in CruiseGraph.surfaces[station_ref].

        The processed surfaces and their labels are also added to CruiseGraph.ctdz in
        columns called "p_at_<station_ref>" and "veronis_<station_ref>".

        Parameters
        ----------
        station_ref : int or str
            The identifier for the station to use as the reference.
        cutoff : float
            The cutoff for outlier removal in dbar, by default 5.  Set to a very high
            number to not remove outliers.
        """
        from neutralocean.surface import omega_surf

        i_ref = self._station_to_index(station_ref)
        surfaces = np.full((self.stations.shape[0], self.stations.npts.max()), np.nan)
        for i, p_init in enumerate(
            self.ctdz[self.ctdz.station == station_ref].pressure
        ):
            # For neutralocean v2.1.3:
            surfaces[:, i] = omega_surf(
                *self.stp,
                self.grid,
                pin_cast=i_ref,
                pin_p=p_init,
                eos="jmdfwg06",
                output=False
            )[2]
            # # For neutralocean v2.2.0:
            # surfaces[:, i] = omega_surf(
            #     *self.stp,
            #     self.grid,
            #     pin_cast=i_ref,
            #     p_init=p_init,
            #     eos="jmdfwg06",
            #     output=False
            # )[2]
        self.surfaces_raw[station_ref] = surfaces
        self._smooth_and_interpolate(station_ref, cutoff=cutoff)
        self._get_ref_surfaces(station_ref)
        self._label_surfaces(station_ref)

    def get_surfaces_all(self, cutoff=5, verbose=True):
        """Compute all omega surfaces using CruiseGraph.get_surfaces() with every
        station used as the reference separately.

        Parameters
        ----------
        cutoff : float
            The cutoff for outlier removal in dbar, by default 5.  Set to a very high
            number to not remove outliers.
        verbose : bool
            Whether to print progress, by default True.
        """
        for i, station_ref in enumerate(self.stations.index):
            if verbose:
                print(
                    "Getting surfaces wrt. station {} ({} of {})...".format(
                        station_ref, i + 1, self.stations.index.size
                    )
                )
            self.get_surfaces(station_ref, cutoff=cutoff)
        if verbose:
            print("Got all surfaces!")

    def _remove_outliers(self, station_ref, cutoff=5):
        surfaces_raw_clean = self.surfaces_raw[station_ref].copy()
        for i, station in enumerate(surfaces_raw_clean):
            for j in range(1, len(station) - 1):
                if (
                    (np.abs(station[j] - station[j - 1]) >= cutoff)
                    & (np.abs(station[j] - station[j + 1]) >= cutoff)
                    & (np.abs(station[j - 1] - station[j + 1]) <= cutoff)
                ):
                    surfaces_raw_clean[i, j] = (station[j - 1] + station[j + 1]) / 2
        return surfaces_raw_clean

    def _smooth_and_interpolate(self, station_ref, cutoff=5):
        surfaces_raw_clean = self._remove_outliers(station_ref, cutoff=cutoff)
        self.surfaces[station_ref] = np.full(surfaces_raw_clean.shape, np.nan)
        pressure_ref = self.get_pressure_ref(station_ref)
        for i in range(len(self.stations.index)):
            level = plot.smooth_whittaker(surfaces_raw_clean[i])
            # Remove NaNs for interpolation and interpolate
            L = ~np.isnan(level)
            interp = interpolate.PchipInterpolator(
                pressure_ref[L], level[L], extrapolate=False
            )
            self.surfaces[station_ref][i] = interp(pressure_ref)

    def _get_ref_surfaces(self, station_ref):
        pressure_ref = self.get_pressure_ref(station_ref)
        self.ctdz["p_at_{}".format(station_ref)] = np.nan
        for i, s in enumerate(self.stations.index):
            L = ~np.isnan(self.surfaces[station_ref][i])
            interp = interpolate.PchipInterpolator(
                self.surfaces[station_ref][i][L], pressure_ref[L], extrapolate=False
            )
            S = self.ctdz.station == s
            self.ctdz.loc[S, "p_at_{}".format(station_ref)] = interp(
                self.ctdz.loc[S, "pressure"]
            )

    def _label_surfaces(self, station_ref):
        ref_veronis = self.ctdz[self.ctdz.station == station_ref][
            ["pressure", "veronis_here"]
        ]
        interp_veronis = interpolate.PchipInterpolator(
            ref_veronis.pressure, ref_veronis.veronis_here
        )
        self.ctdz["veronis_{}".format(station_ref)] = interp_veronis(
            self.ctdz["p_at_{}".format(station_ref)]
        )

    def _station_to_index(self, station):
        return (self.stations.index == station).nonzero()[0][0]

    def _index_to_station(self, index):
        return self.stations.index[index]

    def get_pressure_ref(self, station_ref):
        i_ref = self._station_to_index(station_ref)
        pressure_ref = self.stp[2][i_ref]
        return pressure_ref

    def _get_pos(self, crs):
        self.crs = crs
        self.pos = {
            i: self.crs.transform_points(
                ccrs.PlateCarree(), row.longitude, row.latitude
            )[0][:2]
            for i, row in self.stations.iterrows()
        }

    def plot_surface_map(
        self,
        station_ref,
        p_surface,
        crs=None,
        cutoff=5,
        dpi=300,
        extent=None,
        pad_latitude=0.2,
        pad_longitude=0.2,
    ):
        """Draw the adjacency graph and pressure on a given surface.

        Parameters
        ----------
        station_ref : int or str
            The identifier for the station to use as the reference.  The method
            CruiseGraph.get_surfaces(station_ref) is run if it hasn't already been.
        p_surface : float
            What pressure at the reference station to plot the surface for.  Does not
            need to match exactly, the closest one will be used.
        crs : Cartopy crs
            A coordinate reference system from cartopy.crs.  If not provided, then
            cartopy.crs.PlateCarree() is used.
        cutoff : float, optional
            The cutoff for outlier removal in dbar, by default 5.  Set to a very high
            number to not remove outliers.  Only used if the method
            CruiseGraph.get_surfaces(station_ref) needs to be run.
        dpi : int, optional
            Figure resolution, by default 300.
        extent : array-like, optional
            Map extent, as [lon_min, lon_max, lat_min, lat_max].  If not provided,
            the extent of the stations plus padding is used (see pad_latitude and
            pad_longitude).
        pad_latitude : float, optional
            Fraction by which to pad the latitude axis if extent not provided, by
            default 0.2.
        pad_longitude : float, optional
            Fraction by which to pad the longitude axis if extent not provided, by
            default 0.2.
        """
        if station_ref not in self.surfaces:
            self.get_surfaces(station_ref, cutoff=cutoff)
        if extent is None:
            lat_min = self.stations.latitude.min()
            lat_max = self.stations.latitude.max()
            lon_min = self.stations.longitude.min()
            lon_max = self.stations.longitude.max()
            lat_diff = lat_max - lat_min
            lon_diff = lon_max - lon_min
            extent = [
                lon_min - lon_diff * pad_longitude,
                lon_max + lon_diff * pad_longitude,
                lat_min - lat_diff * pad_latitude,
                lat_max + lat_diff * pad_latitude,
            ]
        pressure_ref = self.get_pressure_ref(station_ref)
        i = np.argmin(np.abs(pressure_ref - p_surface))
        if i >= self.stations.loc[station_ref].npts:
            i = self.stations.loc[station_ref].npts - 1
        i = int(i)
        if crs is None:
            crs = ccrs.PlateCarree()
        self._get_pos(crs)
        fig, ax = plt.subplots(dpi=dpi, subplot_kw={"projection": crs})
        nx.draw_networkx(
            self,
            pos=self.pos,
            ax=ax,
            # Node marker properties
            node_size=0,
            linewidths=0,
            # Edge properties
            width=0.5,
            # Label properties
            font_size=8,
            with_labels=False,
        )
        sc = ax.scatter(
            "longitude",
            "latitude",
            data=self.stations,
            transform=ccrs.PlateCarree(),
            s=20,
            c=self.surfaces[station_ref][:, i],
            zorder=10,
            cmap="viridis_r",
            vmin=self.ctdz.loc[self.ctdz.station == station_ref].pressure.min(),
            vmax=self.ctdz.loc[self.ctdz.station == station_ref].pressure.max(),
        )
        plt.colorbar(sc, label="Neutral surface pressure / dbar")
        ax.scatter(
            "longitude",
            "latitude",
            data=self.stations.loc[station_ref],
            s=60,
            c="none",
            edgecolor="xkcd:strawberry",
            zorder=11,
            transform=ccrs.PlateCarree(),
        )
        ax.set_extent(extent, crs=ccrs.Geodetic())
        ax.set_title(
            "Ref. pressure = {:.1f} dbar".format(
                self.ctdz.loc[self.ctdz.station == station_ref].pressure.iloc[i]
            )
        )
        ax.add_feature(
            cfeature.NaturalEarthFeature(
                "physical", "land", "10m", edgecolor="none", facecolor="xkcd:dark"
            )
        )
        plot.add_credit(ax)
        fig.tight_layout()
        return fig, ax
