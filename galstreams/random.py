
import astropy.coordinates as coord
import astropy.units as u
import numpy as np


__all__ = ['get_uniform_spherical_angles', 'get_uniform_sphere']


@u.quantity_input(lon_lim=u.deg, lat_lim=u.deg)
def get_uniform_spherical_angles(size=1,
                                 lon_lim=[0., 360]*u.deg,
                                 lat_lim=[-90., 90]*u.deg):
    """Generate uniform random positions on the sphere

    Parameters
    ----------
    size : int
        The number of points to generate.
    lon_lim : `~astropy.units.Quantity`
        The longitude limits to generate as an astropy Angle object or Quantity
        with angular units.
    lat_lim : `~astropy.units.Quantity`
        The latitude limits to generate as an astropy Angle object or Quantity
        with angular units.

    Returns
    -------
    representation : `~astropy.coordinates.UnitSphericalRepresentation`
        An astropy unit spherical representation object containing the random
        spherical positions.
    """

    lon = np.random.uniform(lon_lim[0].value,
                            lon_lim[1].value,
                            size) * lon_lim.unit

    K = np.sin(lat_lim[1]) - np.sin(lat_lim[0])
    arg = K * np.random.uniform(size=size) + np.sin(lat_lim[0])
    lat = np.arcsin(arg)

    return coord.UnitSphericalRepresentation(lon, lat)

@u.quantity_input(lon_lim=u.deg, lat_lim=u.deg, dist_lim=[u.one, u.pc])
def get_uniform_sphere(size,
                       lon_lim=[0., 360]*u.deg,
                       lat_lim=[-90., 90]*u.deg,
                       dist_lim=[0, 1.]*u.one):
    """Generate uniform random positions inside a spherical volume.

    i.e. this can be used to generate points uniformly distributed through a
    spherical annulus by specifying the distance limits.

    Parameters
    ----------
    size : int
        The number of points to generate.
    lon_lim : `~astropy.units.Quantity`
        The longitude limits to generate as an astropy Angle object or Quantity
        with angular units.
    lat_lim : `~astropy.units.Quantity`
        The latitude limits to generate as an astropy Angle object or Quantity
        with angular units.
    dist_lim : `~astropy.units.Quantity`
        The distance limits to generate as an astropy Quantity, either
        dimensionless or with length units.

    Returns
    -------
    representation : `~astropy.coordinates.SphericalRepresentation`
        An astropy spherical representation object containing the random
        spherical positions.
    """

    # R distributed as R^2
    r = np.cbrt(np.random.uniform(dist_lim[0].value**3,
                                  dist_lim[1].value**3,
                                  size=size)) * dist_lim.unit

    rep = get_uniform_spherical_angles(size=size,
                                       lon_lim=lon_lim,
                                       lat_lim=lat_lim)

    return coord.SphericalRepresentation(lon=rep.lon,
                                         lat=rep.lat,
                                         distance=r)
