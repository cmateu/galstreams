import numpy as np
import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
import os, os.path
import sys
import astropy.table
import astropy.coordinates as ac
import astropy.units as u
import gala
##Yes, I know this is not particularly pretty, but its legible for me.
#try:
# import gala
# HAS_GALA = True
#except ImportError:
#    HAS_GALA = False
#if not HAS_GALA:
#    raise ImportError(
#            "\n"
#            "galstreams dynamics requires the 'gala' package.\n\n"
#            "To fix this, install it with:\n\n"
#            "    pip install galstreams[dyn]\n\n"
#            "or install directly:\n\n"
#            "    pip install gala\n\n"
#            "If you are using conda, prefer:\n\n"
#            "    conda install -c conda-forge gala\n"
#    )
import gala.coordinates as gc
#import gala.dynamics as gd
from packaging import version

from .track6d import Track6D
from .mwstreams import MWStreams


#---------------------------------
def get_random_spherical_angles(n,az=[0.,2*np.pi],lat=[-np.pi/2.,np.pi/2],degree=False):

    if degree: f=np.pi/180.
    else: f=1.

    rao,raf=az[0]*f,az[1]*f
    deco,decf=lat[0]*f,lat[1]*f
    #Alpha uniformly distributed
    alpha_s=(raf-rao)*np.random.random(n) + rao
    #Delta distributed as cos(delta)
    K=np.sin(decf)-np.sin(deco)
    x=np.random.random(n)
    delta_s=np.arcsin(K*x+np.sin(deco))

    if degree: alpha_s,delta_s=alpha_s*180./np.pi,delta_s*180./np.pi

    return (alpha_s,delta_s)


def get_random_spherical_coords(n,rad=[0.,1.],az=[0.,2*np.pi],lat=[-np.pi/2.,np.pi/2],degree=False):

    #R distributed as R^2
    ro,rf=rad
    x=np.random.random(n)
    C=(1./3.)*(rf**3-ro**3)
    R_s=(3*C*x+ro**3)**(1./3.)

    phi_s,theta_s=get_random_spherical_angles(n,degree=degree,az=az,lat=lat)

    return (R_s,phi_s,theta_s)

def get_avg_vec(phis,thetas,degree=True,lon0=0.):

    X,Y,Z=bovyc.lbd_to_XYZ(phis,thetas,np.ones_like(phis),degree=degree).T

    #Vector sum
    Xsum=X.sum()
    Ysum=Y.sum()
    Zsum=Z.sum()

    #Normalize (not necessary, but nicer)
    norm=np.sqrt(Xsum**2+Ysum**2+Zsum**2)
    Xsum,Ysum,Zsum=Xsum/norm,Ysum/norm,Zsum/norm

    #Back to spherical
    phisum,thetasum,Rgal=bovyc.XYZ_to_lbd(Xsum,Ysum,Zsum,degree=degree)

    return(phisum,thetasum)

def skycoord_to_string(skycoord):

    """ Convert a one-dimenstional list of SkyCoord to string for Gaia's query format (from DataCarpentry)"""
    corners_list_str = skycoord.to_string()
    corners_single_str = ' '.join(corners_list_str)
    return corners_single_str.replace(' ', ', ')

def get_adql_query_from_polygon(skycoo, base_query=None):

    """ Print part of ADQL that selects points inside input polygon given by SkyCoord object

        Parameters:

        base_query : the base ADQL code for your query to have *before* the polygon selection part

    """

    #if dn<1: print ('Invalid N, N=3 is the minimum allowed number of vertices for polygon')

    skycoord_poly = skycoo.transform_to(ac.ICRS)

    sky_point_list = skycoord_to_string(skycoord_poly)
    polygon_query_base = """{base_query}
    1 = CONTAINS(POINT(ra, dec),
                 POLYGON({sky_point_list}))
    """

    #Make sure to warn user if the polygon is too large for the Gaia Archive query to take it
    length =  np.sum(skycoo[0:-1].separation(skycoo[1:]))
    if length > 22.*u.deg: print('WARNING: Gaia Archive ADQL queries do not support polygons longer than 23deg')

    return polygon_query_base.format(base_query=base_query,sky_point_list=sky_point_list)


def plot_5D_tracks_subplots_row(coo , frame, axs=None, name=None, plot_flag='111', scat_kwds=None, show_ylabels=True,
                                show_xlabel=True, show_legend=False, InfoFlags='1111'):


    fr = frame

    #Get representation names for selected frame and flip around
    l = coo.transform_to(fr).representation_component_names
    n = dict((v,k) for k,v in l.items())
    pm1_name = 'pm_{lon}_cos{lat}'.format(lon=n['lon'],lat=n['lat'])
    pm2_name = 'pm_{lat}'.format(lat=n['lat'])

    if axs is None:
        fig, axs = plt.subplots(1,4, figsize=(12,3))
        plt.tight_layout(pad=1.1, w_pad=1.4)


    ax = axs[0]
    axd = axs[1]
    axpm1 = axs[2]
    axpm2 = axs[3]

    if scat_kwds is None: scat_kwds=dict(marker='.', alpha=0.5)

    ax.scatter( getattr(coo.transform_to(fr), n['lon']).value, getattr(coo.transform_to(fr), n['lat']).value,
               **scat_kwds)

    if plot_flag[1]=='1' and InfoFlags[1]!='0': #if distance info not available (set to 0), don't plot it
        axd.scatter( getattr(coo.transform_to(fr), n['lon']).value, getattr(coo.transform_to(fr), n['distance']).value,
                    **scat_kwds)
    if plot_flag[2]=='1' and InfoFlags[2]!='0': #if pm info not available (set to 0), don't plot it
        axpm1.scatter( getattr(coo.transform_to(fr), n['lon']).value,
                    getattr(coo.transform_to(fr), pm1_name).value,
                    **scat_kwds)
        axpm2.scatter( getattr(coo.transform_to(fr), n['lon']).value,
                    getattr(coo.transform_to(fr), pm2_name).value,
                    **scat_kwds)

    if show_legend: ax.legend()#, bbox_to_anchor=(1.02,0.95))

    if show_xlabel:
        for ax_i in [ax,axd,axpm1,axpm2]:
            ax_i.set_xlabel("{lon} ({unit})".format(lon=n['lon'], unit=getattr(coo.transform_to(fr), n['lon']).unit))

    if show_ylabels:
        ax.set_ylabel("{y} ({unit})".format(y=n['lat'], unit=getattr(coo.transform_to(fr), n['lat']).unit))
        axpm1.set_ylabel("{y} ({unit})".format(y=pm1_name, unit=getattr(coo.transform_to(fr), pm1_name).unit))
        axpm2.set_ylabel("{y} ({unit})".format(y=pm2_name, unit=getattr(coo.transform_to(fr), pm2_name).unit))
        axd.set_ylabel("D (kpc)")


def get_mask_in_poly_footprint(poly_sc, coo, stream_frame):

     ''' Test whether points in input SkyCoords object are inside polygon footprint. 

         Parameters
         ==========

         poly_sc : astropy.coordinates.SkyCoord object with polygon vertices 
         coo : astropy.coordinates.SkyCoord object

         Returns
         =======

         mask : boolean mask array, same number of elements as coo 
     '''

     #Create poly-path object
     verts = np.array([poly_sc.transform_to(stream_frame).phi1, poly_sc.transform_to(stream_frame).phi2]).T
     poly = mpl.path.Path(verts)

     #The polygon test has to be done in phi1/phi2 (otherwise there's no guarantee of continuity for the polygon)
     coo_in_str_fr = coo.transform_to(stream_frame)
     _points = np.stack((coo_in_str_fr.phi1, coo_in_str_fr.phi2)).T

     return poly.contains_points(_points)

def _get_gala_dynamics():
    try:
        import gala.dynamics as gd
        return gd
    except Exception as err:
        raise ImportError(
            "\n"
            "galstreams requires a working gala installation.\n\n"
            "The current gala install is broken (likely GSL / compiled dependency issue).\n\n"
            "Fix options:\n\n"
            "  pip install --force-reinstall gala\n"
            "or (recommended):\n"
            "  conda install -c conda-forge gala gsl\n"
        ) from err

def compute_angular_momentum_track(track, return_cartesian = False):

   '''  Compute angular momentum for each point in the track.  
        By default it returns the spherical components of the angular momentum in the heliocentric and galactocentric reference
        frames at rest w.r.t. the GSR. Set return_cartesian = True to get cartesian components
 

       	Parameters:
      	=======

	track : SkyCoord object
      
        return_cartesian : If True returns cartesian coordinates. If False, returns spherical coords (astropy format mod, lat, lon)
     
        Returns:
	========

        L : list object with compoments of angular momentum vector. By default returns spherical components modulus, lat, lon

   '''

   gd = _get_gala_dynamics() 

   tr = track.cartesian

   pos = ac.CartesianRepresentation(x = tr.x, y = tr.y, z = tr.z)
   vel = ac.CartesianDifferential(d_x = tr.differentials['s'].d_x, d_y = tr.differentials['s'].d_y, d_z = tr.differentials['s'].d_z)
   psp = gd.PhaseSpacePosition(pos=pos, vel=vel)
   L = psp.angular_momentum()

   if return_cartesian: return L
   else: 
      L_sph = ac.cartesian_to_spherical(x = L[0], y = L[1], z = L[2])
      return L_sph


#-----------------------------------------------------------------------------------------------

def plot_globular_clusters(ax,plot_colorbar=False,scat_kwargs=None,galactic=True):

    #Harris's Globular cluster compilation
    gcfilen=os.path.join(os.path.dirname(os.path.realpath(__file__)),'lib','globular_cluster_params.harris2010.tableI.csv')
    gc_l,gc_b,gc_Rhel=scipy.genfromtxt(gcfilen,missing_values="",usecols=(5-1,6-1,7-1),unpack=True,delimiter=',',dtype=np.float)

    gc_RA,gc_DEC=bovyc.lb_to_radec(gc_l,gc_b,degree=True).T

    #set a few plotting and labelling defaults
    scatter_kwargs=dict(marker='x',s=70.,vmin=0.,zorder=0)
    #but override whichever are user-supplied (doing it this way I ensure using my own defaults and not matplotlib's
    #if user supplies values for some (but not all) keywords
    if scat_kwargs is not None:
        for key in scat_kwargs.keys(): scatter_kwargs[key]=scat_kwargs[key]

    #Plot globular cluster layer
    if galactic: cc=ax.scatter(gc_l,gc_b,c=gc_Rhel,**scatter_kwargs)
    else: cc=ax.scatter(gc_RA,gc_DEC,c=gc_Rhel,**scatter_kwargs)

    if plot_colorbar: plt.colorbar(cc,ax=ax)

#----------------
