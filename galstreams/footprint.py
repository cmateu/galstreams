# Third-party
import astropy.coordinates as co
import astropy.units as u
import numpy as np
from spherical_geometry.polygon import SingleSphericalPolygon

from .config import config
from .random import get_uniform_spherical_angles

class StreamFootprint:

    def __init__(self, name, poly, frame=None):
        """Represents the sky footprint of a given stellar stream.

        Parameters
        ----------
        name : str
            The stream name.
        polygon : `~astropy.units.Quantity`, `~spherical_geometry.SphericalPolygon`
            Either an astropy quantity object (with angular units) specifying
            the vertices of a spherical polygon, or a spherical polygon
            instance.
        frame : `astropy.coordinates.BaseCoordinateFrame` subclass instance (optional)
            If not specified, the default is taken from the config settings.
        """

        self.name = str(name)

        # Validate input polygon. Can be a set of vertices, or an object
        if not isinstance(poly, SingleSphericalPolygon):
            poly = u.Quantity(poly)
            if poly.ndim != 2 or poly.shape[1] != 2:
                raise ValueError("Invalid shape for polygon vertices {}. "
                                 "Expected something like (N, 2) where N is "
                                 "the number of vertices.".format(poly.shape))
            poly = SingleSphericalPolygon.from_lonlat(
                poly[:, 0].to_value(u.deg), poly[:, 1].to_value(u.deg),
                degrees=True)

        self.poly = poly

        # If no frame is specified, assume that the input polygon is in a frame
        # specified in the configuration.
        if frame is None:
            frame = config['default_frame']

        elif hasattr(frame, 'frame'): # if a SkyCoord instance
            frame = frame.frame

        if not isinstance(frame, co.BaseCoordinateFrame):
            raise ValueError("Input coordinate frame must be an astropy "
                             "coordinates frame subclass *instance*, not a "
                             "'{}'".format(frame.__class__.__name__))
        self.frame = frame

        # Store an internal ICRS representation of the footprint
        lon, lat = self.poly.to_lonlat()
        self._poly_c = co.SkyCoord(lon * u.deg, lat * u.deg,
                                   frame=self.frame)

    def from_corners(self, corners):
        """ TODO: """
        pass

    def from_galactocentric(self, name, poly_c, frame=None):
        """ TODO: """
        pass

    def get_polygon(self, frame=None):
        """Get a ``SingleSphericalPolygon`` object in the specified coordinate
        frame to represent the footprint.

        Parameters
        ----------
        frame : `astropy.coordinates.BaseCoordinateFrame` subclass instance (optional)
            If not specified, the default is taken from the config settings.

        Returns
        -------
        poly : ``spherical_geometry.SingleSphericalPolygon``
            The polygon object.

        """
        if frame is None:
            frame = config['default_frame']

        poly_c = self._poly_c.transform_to(frame)
        poly = SingleSphericalPolygon.from_lonlat(
            poly_c.spherical.lon.degree,
            poly_c.spherical.lat.degree)

        return poly

    def get_sky_center(self, frame=None):
        """Estimate the central point of the polygon

        This estimate will become worse when the number of vertices along the
        polygon is small.

        Parameters
        ----------
        frame : `astropy.coordinates.BaseCoordinateFrame` subclass instance (optional)
            If not specified, the default is taken from the config settings.

        Returns
        -------
        midpt : `~astropy.coordinates.SkyCoord`
            The coordinates of the polygon center.

        """
        poly = self.get_polygon(frame)

        midpt = co.CartesianRepresentation(*np.mean(poly.points, axis=0)*u.kpc)\
            .represent_as(co.UnitSphericalRepresentation)

        return co.SkyCoord(midpt, frame=frame)

    def get_random_points(self, size=1024, frame=None, random_state=None,
                          wrap_angle=360*u.deg, max_niter=128):
        """Generate uniform random points throughout the stream footprint

        Parameters
        ----------
        size : int (optional)
            The number of points to generate.
        frame : `astropy.coordinates.BaseCoordinateFrame` subclass instance (optional)
            The coordinate frame to generate points in.
        random_state : `numpy.random.RandomState` (optional)
            A numpy random state object used to control the random number
            generator and seed.
        wrap_angle : `~astropy.coordinates.Angle` (optional)
            The angle to wrap the longitude coordinate at. Specify this if you
            are having issues with the polygon wrapping around 360 deg.
        max_niter : int (optional)
            The maximum number of iterations to use when iteratively generating
            points within the footprint.

        Returns
        -------
        c : `~astropy.coordinates.Skycoord`
            The coordinates of random points within the sky footprint.

        """
        if frame is None:
            frame = config['default_frame']

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
            raise ValueError("Exceeded maximum number of iterations when "
                             "generating points. Try increasing `max_niter` to "
                             "something like max_niter=1024 (but note this "
                             "might take a while to run!).")

        rep = co.concatenate_representations(reps)
        return co.SkyCoord(rep[:size], frame=frame)

    def contains_coords(self, c):
        """Check whether the input coordinate or coordinates are inside of the
        stream footprint.

        Parameters
        ----------
        c : `~astropy.coordinates.SkyCoord`
            Either a scalar coordinate or a set of coordinates to check.

        Returns
        -------
        in_footprint : bool, `numpy.ndarray`
            Either a boolean (for scalar input) or a boolean array (for
            multiple coordinates) specifying which of the input coordinates are
            inside of the stream footprint.
        """
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

    def plot_polygon(self, frame=None, ax=None,
                     wrap_angle=360*u.deg, autolim=True, **kwargs):
        """Plot the footprint polygon

        Parameters
        ----------
        frame : `~astropy.coordinates.BaseCoordinateFrame` subclass instance (optional)
            The coordinate frame to generate points in.
        ax : `~matplotlib.axes.Axes` (optional)
            The matplotlib axes object to plot on. If not specified, this will
            plot on the current axes object.
        wrap_angle : `~astropy.coordinates.Angle` (optional)
            The angle to wrap the longitude coordinate at. Specify this if you
            are having issues with the polygon wrapping around 360 deg.
        autolim : bool (optional)
            Determine and set the axes limits automatically to show the whole
            polygon.
        **kwargs
            All other keyword arguments are passed to
            `~matplotlib.patches.Polygon`.

        Returns
        -------
        fig : `~matplotlib.figure.Figure`
        ax : `~matplotlib.axes.Axes`
        
        """
        import matplotlib as mpl

        poly = self.get_polygon(frame=frame)
        verts = np.stack(poly.to_lonlat()).T
        verts[:, 0] = co.Angle(verts[:, 0]*u.deg).wrap_at(wrap_angle).degree

        if ax is None:
            import matplotlib.pyplot as plt
            ax = plt.gca()

        fig = ax.figure

        # TODO: Detect if the polygon is mangled by wrapping issues by looking
        # for discontinuous jumps of, e.g., >120 degrees in successive points
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
