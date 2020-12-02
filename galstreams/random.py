
import astropy.coordinates as coord
import astropy.units as u
import numpy as np


__all__ = ['get_uniform_spherical_angles', 'get_uniform_sphere']


@u.quantity_input(lon_lim=u.deg, lat_lim=u.deg)
def get_uniform_spherical_angles(size=1,
                                 lon_lim=[0., 360]*u.deg,
                                 lat_lim=[-90., 90]*u.deg,
                                 random_state=None):
    """Generate uniform random positions on the sphere

    Parameters
    ----------
    size : int
        The number of points to generate.
    lon_lim : `~astropy.units.Quantity` (optional)
        The longitude limits to generate as an astropy Angle object or Quantity
        with angular units.
    lat_lim : `~astropy.units.Quantity` (optional)
        The latitude limits to generate as an astropy Angle object or Quantity
        with angular units.
    random_state : `numpy.random.RandomState` (optional)
        A numpy random state object used to control the random number generator
        and seed.

    Returns
    -------
    representation : `~astropy.coordinates.UnitSphericalRepresentation`
        An astropy unit spherical representation object containing the random
        spherical positions.
    """
    if random_state is None:
        random_state = np.random

    lon = np.random.uniform(lon_lim[0].value,
                            lon_lim[1].value,
                            size) * lon_lim.unit

    K = np.sin(lat_lim[1]) - np.sin(lat_lim[0])
    arg = K * random_state.uniform(size=size) + np.sin(lat_lim[0])
    lat = np.arcsin(arg)

    return coord.UnitSphericalRepresentation(lon, lat)


@u.quantity_input(lon_lim=u.deg, lat_lim=u.deg, dist_lim=[u.one, u.pc])
def get_uniform_sphere(size,
                       lon_lim=[0., 360]*u.deg,
                       lat_lim=[-90., 90]*u.deg,
                       dist_lim=[0, 1.]*u.one,
                       random_state=None):
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
    random_state : `numpy.random.RandomState` (optional)
        A numpy random state object used to control the random number generator
        and seed.

    Returns
    -------
    representation : `~astropy.coordinates.SphericalRepresentation`
        An astropy spherical representation object containing the random
        spherical positions.
    """
    if random_state is None:
        random_state = np.random

    rep = get_uniform_spherical_angles(size=size,
                                       lon_lim=lon_lim,
                                       lat_lim=lat_lim,
                                       random_state=random_state)

    # R distributed as R^2
    r = np.cbrt(random_state.uniform(dist_lim[0].value**3,
                                     dist_lim[1].value**3,
                                     size=size)) * dist_lim.unit

    return coord.SphericalRepresentation(lon=rep.lon,
                                         lat=rep.lat,
                                         distance=r)
