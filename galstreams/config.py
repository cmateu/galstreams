# Third-party
from astropy.coordinates import Galactocentric, CartesianDifferential
import astropy.units as u

config = dict()
config['galcen_distance'] = 8.1 * u.kpc
config['galcen_v_sun'] = [11.1, 232.24, 7.25] * u.km/u.s
config['z_sun'] = 27 * u.pc

def get_galactocentric_frame():
    return Galactocentric(
        galcen_distance=config['galcen_distance'],
        galcen_vsun=CartesianDifferential(config['galcen_v_sun']))


# Default frame for showing:
config['default_frame'] = 'icrs'
