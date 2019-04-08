# Third-party
import astropy.coordinates as co
import astropy.units as u
import numpy as np
from spherical_geometry.polygon import SingleSphericalPolygon

from .random import get_uniform_spherical_angles

class StreamFootprint:

    @u.quantity_input(distance=u.kpc)
    def __init__(self, name, poly, frame='icrs', distance=None):
        """TODO: describe

        Parameters
        ----------
        name : str
        polygon : `spherical_geometry.SphericalPolygon`
        frame : `astropy.coordinates.BaseCoordinateFrame` subclass instance (optional)
            Default is ICRS.
        """

        self.name = str(name)
        self.sname = self.name[:8] # TODO: ???

        # TODO: validate polygon. Can be list of vertices, or object
        self.poly = poly

        if isinstance(frame, str):
            frame = co.frame_transform_graph.lookup_name(frame)
        self.frame = frame

        if isinstance(self.frame, co.Galactocentric) and distance is None:
            raise ValueError("TODO")

        # Store an internal ICRS representation of the footprint
        lon, lat = self.poly.to_lonlat()
        if distance is not None:
            poly_c = co.SkyCoord(lon * u.deg, lat * u.deg,
                                 distance=distance,
                                 frame=self.frame)
        else:
            poly_c = co.SkyCoord(lon * u.deg, lat * u.deg,
                                 frame=self.frame)
        self._poly_c = poly_c

        # Store ICRS poly in cache by default (this triggers the caching):
        self.get_polygon(frame=co.ICRS)

    def from_corners(self, corners):
        """ TODO: """
        pass

    def get_polygon(self, frame='icrs'):
        """TODO:

        Parameters
        ----------
        frame : str, `~astropy.coordinates.BaseCoordinateFrame`
        """
        if not isinstance(frame, str):
            frame = frame.name

        poly_c = self._poly_c.transform_to(frame)
        poly = SingleSphericalPolygon.from_lonlat(
            poly_c.spherical.lon.degree,
            poly_c.spherical.lat.degree)

        return poly

    def get_midpoint(self, frame='icrs'):
        poly = self.get_polygon(frame)

        midpt = co.CartesianRepresentation(*np.mean(poly.points, axis=0)*u.kpc)\
            .represent_as(co.UnitSphericalRepresentation)

        return co.SkyCoord(midpt, frame=frame)

    def get_random_points(self, size=1024, frame='icrs', random_state=None,
                          wrap_angle=360*u.deg, max_niter=128):
        poly = self.get_polygon(frame=frame)
        verts = np.stack(poly.to_lonlat()).T
        verts[:, 0] = co.Angle(verts[:, 0]*u.deg).wrap_at(wrap_angle).degree

        lon_lim = [verts[:, 0].min(), verts[:, 0].max()] * u.deg
        lat_lim = [verts[:, 1].min(), verts[:, 1].max()] * u.deg

        random_size = 4 * size # MAGIC NUMBER
        N = 0
        niter = 0
        reps = []
        while niter < max_niter:
            random_size = 2 * random_size
            rep = get_uniform_spherical_angles(size=random_size,
                                               lon_lim=lon_lim,
                                               lat_lim=lat_lim,
                                               random_state=random_state)
            c = co.SkyCoord(rep, frame=frame)
            mask = self.contains_coord(c)
            reps.append(rep[mask])
            N += mask.sum()

            if N >= size:
                break

        else:
            raise ValueError("TODO increase max_niter")

        rep = co.concatenate_representations(reps)
        return co.SkyCoord(rep[:size], frame=frame)

    def contains_coord(self, c):
        if isinstance(c, co.BaseRepresentation):
            poly = self.get_polygon()
            rep = c.represent_as(co.SphericalRepresentation)
        else:
            c = co.SkyCoord(c)
            poly = self.get_polygon(frame=c.frame)
            rep = c.spherical

        if c.isscalar:
            return poly.contains_lonlat(rep.lon.degree,
                                        rep.lat.degree)

        else:
            mask = [poly.contains_lonlat(x, y)
                    for x, y in zip(rep.lon.degree,
                                    rep.lat.degree)]
            return np.array(mask)

    def plot_polygon(self, frame='icrs', ax=None,
                     wrap_angle=360*u.deg, autolim=True, **kwargs):
        """TODO"""
        import matplotlib as mpl

        poly = self.get_polygon(frame=frame)
        verts = np.stack(poly.to_lonlat()).T
        verts[:, 0] = co.Angle(verts[:, 0]*u.deg).wrap_at(wrap_angle).degree

        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        # Detect if the polygon is mangled by wrapping issues by looking
        # for discontinuous jumps of >120 degrees in successive points
        # - see APW's notebook Galstreams.ipynb
        if np.any(np.abs(np.diff(verts[:, 0])) > 120):
            # TODO: for now, raise a warning
            pass

        poly_patch = mpl.patches.Polygon(verts, **kwargs)
        ax.add_patch(poly_patch)

        if autolim:
            h = 0.05 # make limits 5% larger than span
            xlim = (verts[:, 0].min(), verts[:, 0].max())
            span = xlim[1]-xlim[0]
            xlim = (xlim[1] + h*span, xlim[0] - h*span)

            ylim = (verts[:, 1].min(), verts[:, 1].max())
            span = ylim[1]-ylim[0]
            ylim = (ylim[0] - h*span, ylim[1] + h*span)

            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

        return fig, ax
