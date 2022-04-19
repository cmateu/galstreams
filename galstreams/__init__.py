import numpy as np
import scipy
import matplotlib as mpl
import pylab as plt
import gcutils 
import os, os.path
import sys
import astropy.table
import astropy.coordinates as ac
import astropy.units as u
import gala.coordinates as gc
import gala.dynamics as gd

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

   tr = track.cartesian

   pos = ac.CartesianRepresentation(x = tr.x, y = tr.y, z = tr.z)
   vel = ac.CartesianDifferential(d_x = tr.differentials['s'].d_x, d_y = tr.differentials['s'].d_y, d_z = tr.differentials['s'].d_z)
   psp = gd.PhaseSpacePosition(pos=pos, vel=vel)
   L = psp.angular_momentum()

   if return_cartesian: return L
   else: 
      L_sph = ac.cartesian_to_spherical(x = L[0], y = L[1], z = L[2])
      return L_sph


#---------MW Streams class--------------------------------------------------------------------------------
class MWStreams(dict):
    
  def __init__(self, verbose=False, implement_Off=False):

    #A MWStreams object is a dictionary in which each entry is a Footprint object, indexed by each stream's name.
    #There's also a mandatory summary entry, which is a Pandas DataFrame with summary attributes for the full library

    #Initialize empty dictionary
    dict.__init__(self)

    #Read in the master logs
    tdir = os.path.dirname(os.path.realpath(__file__))
    #master logs
    filepath = "{path}/{filen}".format(path=tdir+"/lib/",filen='master_log.txt')
    lmaster = astropy.table.Table.read(filepath,format='ascii.commented_header').to_pandas()
    filepath = "{path}/{filen}".format(path=tdir+"/lib/",filen='master_log.discovery_refs.txt')
    lmaster_discovery = astropy.table.Table.read(filepath,format='ascii.commented_header').to_pandas(index='Name')
    filepath = "{path}/{filen}".format(path=tdir+"/lib/",filen='master_log.comments.txt')
    lmaster_comments = astropy.table.Table.read(filepath,format='ascii.commented_header').to_pandas(index='Name')

    lmaster["On"] = lmaster["On"].astype('bool') #this attribute controls whether a given track is included or not

    #SkyCoords objects will be created for each of these dicts after the full library has been created
    attributes = ['ra','dec','distance','pm_ra_cosdec','pm_dec','radial_velocity']
    units = [u.deg, u.deg, u.kpc, u.mas/u.yr, u.mas/u.yr, u.km/u.s]
    end_o_dic = {k: np.array([])*uu  for k,uu in zip(attributes,units) }
    end_f_dic = {k: np.array([])*uu  for k,uu in zip(attributes,units) }
    mid_point_dic = {k: np.array([])*uu  for k,uu in zip(attributes,units) }
    mid_pole_dic  = {k: np.array([])*uu  for k,uu in zip(attributes,units) }
    info_flags = []
    discovery_refs = []
    lengths = np.array([])*u.deg

    print("Initializing galstreams library from master_log... ")
    nid = 1
    for ii in np.arange(lmaster.TrackRefs.size):

       #Create the names of the files containing the knots and summary attributes to initialize each stream
       summary_file = "{tdir}/track.{imp}.{stname}.{tref}.summary.ecsv".format(tdir=tdir+'/tracks', imp=lmaster.Imp[ii], 
          								      stname=lmaster.Name[ii],
          								      tref=lmaster.TrackRefs[ii])
       track_file   = "{tdir}/track.{imp}.{stname}.{tref}.ecsv".format(tdir=tdir+'/tracks', imp=lmaster.Imp[ii], 
          							      stname=lmaster.Name[ii],
          							      tref=lmaster.TrackRefs[ii])

       if verbose: 
          print(f"Initializing Track6D {lmaster.TrackName[ii]} for {lmaster.Name[ii]}...")

       #Do the magic. The track is read and all attributes stored in the summary for all registered stream tracks. 
       #Only the ones "On" are "realized" unless implement_Off == True
       track = Track6D(track_name=lmaster.TrackName[ii], stream_name=lmaster.Name[ii], track_reference=lmaster.TrackRefs[ii], 
                       track_file=track_file, track_discovery_references=lmaster_discovery.loc[lmaster.Name[ii],'DiscoveryRefs'] ,
                       summary_file=summary_file)

#       self[lmaster.TrackName[ii]] = track
       if implement_Off:
           self[lmaster.TrackName[ii]] = track
           self[lmaster.TrackName[ii]].ID = nid
           nid = nid+1
       else:
          if lmaster.On[ii]: 
              self[lmaster.TrackName[ii]] = track
              self[lmaster.TrackName[ii]].ID = nid
              nid = nid+1
          elif verbose: print(f"Skipping Off track {lmaster.TrackName[ii]}...")

       #Store summary attributes
       for k in attributes: end_o_dic[k] = np.append(end_o_dic[k], getattr(track.end_points, k)[0] )
       for k in attributes: end_f_dic[k] = np.append(end_f_dic[k], getattr(track.end_points, k)[1] )
       for k in attributes: mid_point_dic[k] = np.append(mid_point_dic[k], getattr(track.mid_point, k) )
       for k in attributes[:2]: mid_pole_dic[k]  = np.append(mid_pole_dic[k] , getattr(track.mid_pole, k) )

       info_flags.append(track.InfoFlags)
       lengths = np.append(lengths, track.length)
       discovery_refs = np.append(discovery_refs, lmaster_discovery.loc[lmaster.Name[ii],'DiscoveryRefs'] )


    #Add skycoord summary attributes to the library and selected cols to the summary table
    self.end_o = ac.SkyCoord(**end_o_dic)
    self.end_f = ac.SkyCoord(**end_f_dic)
    self.mid_point = ac.SkyCoord(**mid_point_dic)
    self.mid_pole  = ac.SkyCoord(ra=mid_pole_dic["ra"], dec=mid_pole_dic["dec"], frame='icrs')

    #Store master table as an attribute (inherits structure of lmaster dataframe)
    self.summary = lmaster

    #Stream Length
    self.summary["length"] = np.array(lengths)
    #End points
    self.summary["ra_o"] = end_o_dic["ra"]
    self.summary["dec_o"] = end_o_dic["dec"]
    self.summary["distance_o"] = end_o_dic["distance"]
    self.summary["ra_f"] = end_f_dic["ra"]
    self.summary["dec_f"] = end_f_dic["dec"]
    self.summary["distance_f"] = end_f_dic["distance"]
    #Mid point
    self.summary["ra_mid"] = mid_point_dic["ra"]
    self.summary["dec_mid"] = mid_point_dic["dec"]
    self.summary["distance_mid"] = mid_point_dic["distance"]
    #Pole
    self.summary["ra_pole"] = mid_pole_dic["ra"]
    self.summary["dec_pole"] = mid_pole_dic["dec"]
    #Info
    self.summary["InfoFlags"] = np.array(info_flags)
    self.summary["DiscoveryRefs"] = discovery_refs
    #Index by TrackName
    self.summary.index=self.summary.TrackName

    #Create a numeric ID for each track
    self.summary["ID"] = ''
    for ii in self.summary.index:
       if self.summary.loc[ii,'On']:
        self.summary.loc[ii,'ID'] = self[ii].ID
  
    #Store discovery references in summary table
    #for ii in self.summary.index:
    #    self.summary.loc[ii,'DiscoveryRefs'] = self[ii].ref_discovery

  def all_active_track_names(self):
       return self.keys()

  def all_unique_stream_names(self):
       return np.unique(self.summary.Name[self.summary.On])

  def all_track_names(self):
       return np.array(self.summary.index)       

  def get_track_names_for_stream(self, StreamName):
    '''
        Parameters
        ==========

        StreamName : str 
                     Name of the stream (or part of it)
  
        Returns 
        =======

        array : Array with all the TrackNames for which the stream's name matches the input string  

    '''

    all_track_names = []

    for tn in self.summary.index:
      if StreamName in self.summary.loc[tn,'Name']: 
         all_track_names.append(tn)

    return all_track_names

  def plot_stream_compilation(ax, plot_colorbar=True, scat_kwargs=None, use_shortnames=False, 
                                  cb_kwargs=None, verbose=False):

      return

class Track6D:

  def __init__(self, track_name, track_file, summary_file, stream_name=None, stream_shortname=None, track_reference=' ', 
	        track_discovery_references=' ', verbose=True):

      ''' Track6D: A Stellar Stream's Track realization in 6D. See the list of attributes below. 
 	
	Parameters
	==========

	track_name : str
	  Unique identifier for the stream's track realization. Not necesarily identical to stream_name, e.g. if
	  more than one track for the same stream is available
	track_file : str
	  Input ecsv file containing knots to initialize 6D track stream realization
	summary_file : str
	  Input ecsv file containing end point, mid point and pole coordinates (6D, 6D and 3D)

        Attributes:
	===========

	track_name : str
	  Unique identifier for the stream's track realization

	stream_name: str
	  Stream name  

	stream_shortname: str  
	  Stream short name  

	track : astropy.coordinates.SkyCoord Object
	  Contains the track 6D info. By default initialized in icrs frame

        length: angular length measured along the track

        InfoFlags: string - 4-bits indicate available (or assumed) data. 
  	    bit 0: 0 = great circle by construction
  	    bit 1: 0 = no distance track available (only mean or central value reported)
  	    bit 2: 0 = no proper motion data available (only mean or central value reported)
  	    bit 3: 0 = no radial velocity data available (only mean or central value reported)
 
	end_points: 2-element astropy.coordinates.SkyCoord Object with end point coordinates
	  
        mid_point: astropy.coordinates.SkyCoord Object with stream's mid-point coordinates (phi1=0)

        mid_pole: astropy.coordinates.SkyCoord Object heliocentric pole at phi1=0
 
        poly_sc: astropy.coordinates.SkyCoord Object containing vertices for stream's polygon footprint
	
  	mid_pole_gsr: astropy.coordinates.SkyCoord Object. GSR pole at phi1=0

	pole_track_helio: astropy.coordinates.SkyCoord Object heliocentric pole track (galactic coordinates by default)
	
        pole_track_gsr: astropy.coordinates.SkyCoord Object GSR pole track (galactic coordinates by default)
	
        angular_momentum_helio: list object with spherical components (modulus, lat, lon) for the angular momentum of
				each point along the track, computed in a heliocentric frame at rest w.r.t. the GSR
      
        WARNING: angular momentum and pole tracks have length track.size-1  

	'''      


      #Initialize a (new) Footprint6D object

      #First the track name (track's identifier, for popular streams there can be more than one track for a given stream)
      self.track_name = str(track_name)

      #Stream's name
      if stream_name is not None: self.stream_name = str(stream_name)
      else: self.stream_name = self.track_name.copy()

      #Stream's short name
      self.stream_sname = str(stream_shortname)
      if stream_shortname is not None: self.stream_shortname = str(stream_shortname)
      else: self.stream_shortname = self.stream_name[:5]

      #References for the track
      self.ref = track_reference
      self.ref_discovery = track_discovery_references

      #Read-in knots and initialize track
      t = astropy.table.QTable.read(track_file)

      #Store the track in attribute
      self.track = ac.SkyCoord(**t)

      #Now read in the stream's summary file 
      sfile=astropy.table.QTable.read(summary_file)

      #Streams InfoFlags: four-bit flag
      # bit 0: 0 = great circle by construction
      # bit 1: 0 = no distance track available (only mean or central value reported)
      # bit 2: 0 = no proper motion data available (only mean or central value reported)
      # bit 3: 0 = no radial velocity data available (only mean or central value reported)
      self.InfoFlags = str(sfile["InfoFlags"][0]) # All-in-one flag

      #And create the end_points object 
      #two-element SkyCoord obj, one for each end
      end_points = dict()
      atts = [x.replace('end_o.','') for x in sfile.keys() if 'end_o' in x ]
      for att in atts:  #we're effectively looping over skycoords defined for end_o here (ra, dec, ...)
         end_points[att] = np.append(sfile[f'end_o.{att}'],sfile[f'end_f.{att}'])
      self.end_points = ac.SkyCoord(**end_points)

      #Mid-point 
      x = dict()
      atts = [x.replace('mid.','') for x in sfile.keys() if 'mid' in x ]
      for att in atts:  #we're effectively looping over skycoords defined for mid here (ra, dec, ...)
         x[att] = sfile[f'mid.{att}'][0]   #<- make sure to set it up as a scalar. if not, frame conversions get into trouble
      self.mid_point = ac.SkyCoord(**x) 

      #Pole at mid point - The track's (approx) pole at the mid-point. It represents the orbital plane's normal 
      #at the midpoint. If the track is not a great circle as seen from the sun this may change significantly along the track
      x = dict()
      atts = [x.replace('pole.','') for x in sfile.keys() if 'pole' in x ]
      for att in atts:  #we're effectively looping over skycoords defined for pole here (ra, dec, ...)
         x[att] = sfile[f'pole.{att}'][0]
      #Make sure to set the pole's distance attribute to 1 (zero causes problems, when transforming to stream frame coords)
      x["distance"] = 1.*u.kpc   #it shouldn't matter, but if it's zero it does crazy things
      self.mid_pole = ac.SkyCoord(**x)

      #Set up stream's coordinate frame
      self.stream_frame = gc.GreatCircleICRSFrame(pole=self.mid_pole, ra0=self.mid_point.icrs.ra)

      #Compute and store polygon vertices
      self.poly_sc = self.create_sky_polygon_footprint_from_track(width=1*u.deg)

      #Compute and store angular momentum track
      self.angular_momentum_helio = self.get_helio_angular_momentum_track()

      #Compute and store heliocentric pole track
      self.pole_track_helio, self.pole_track_gsr = self.get_pole_tracks()

      #Also store the mid_pole_gsr (poles do not transform as normal coord objects, so this needs to be computed at the GSR)
      #I use this shortcut. The midpoint is located at (helio-centric) phi1=0, so we can retrieve its pole in the gsr track
      mask = np.argmin(np.abs(self.track.transform_to(self.stream_frame).phi1.deg)) #Find central point (closest to phi1=0)
      self.mid_pole_gsr = self.pole_track_gsr[mask]

      #Compute and store the full length along the track
      self.length = np.sum(self.track[0:-1].separation(self.track[1:]))


  def get_helio_angular_momentum_track(self, return_cartesian = False ):

     '''Compute angular momentum for each point in the track.  

        By default it returns the spherical components of the angular momentum in the heliocentric and galactocentric reference
        frames at rest w.r.t. the GSR. If return_cartesian = True it will return cartesian components
     '''

     st_s = self.track.galactic
     #If I wanted the GSR ang momentum: (for now, it doesn't make sense to provide this, it will be misleading as there
     #are too few streams that have radial velocity data.
     #tr = st_s.transform_to(gsr)

     #Heliocentric, at rest w.r.t. GSR ( r_helio x v_gsr ). In this frame the radial velocity component of the stream (and the Sun)
     #does not contribute to the angular momentum. 
     tr = gc.reflex_correct(st_s)

     #Force it to zero if the track doesn't have pm data
     if self.InfoFlags[2]=='0': 
        zz = np.zeros(st_s.b.size) 
        if return_cartesian: 
           L = (zz*u.kpc*u.km/u.s, zz*u.kpc*u.km/u.s, zz*u.kpc*u.km/u.s)
        else:  
           L = (zz*u.kpc*u.km/u.s, zz*u.deg, zz*u.deg)
     else:
        L = compute_angular_momentum_track(tr, return_cartesian = return_cartesian)

     if return_cartesian: return L 
     else:
        return (L[0], L[1].to(u.deg), L[2].to(u.deg) )  #Force lat,lon to be in deg because leaving them in rad is asking for trouble


  def get_pole_tracks(self, use_gsr_default=True):

     ''' TODO '''

     if use_gsr_default: _ = ac.galactocentric_frame_defaults.set('latest')

     #The pole_from_endpoints only works with SkyCoords objs that have no differentials (i.e. no pm/vrad)
     ep1 = ac.SkyCoord(ra=self.track.ra[:-1], dec=self.track.dec[:-1], distance=self.track.distance[:-1], frame='icrs') 
     ep2 = ac.SkyCoord(ra=self.track.ra[1:],  dec=self.track.dec[1:],  distance=self.track.distance[1:], frame='icrs') 

     #That's it. Really, that's it. I love gala. Thanks APW.
     pole_track_helio = gc.pole_from_endpoints(ep1,ep2)
     pole_track_helio = ac.SkyCoord(ra=pole_track_helio.ra, dec=pole_track_helio.dec, frame='icrs')
     #Recast as SkyCoord object (output from prev is ICRS obj, this way coord transf are easier)
     #and make the pole tracks stay in only one hemisphere (we don't have the sense of rotation info in these anyway)
     #the flipping is done in galactic coords (doesn't make sense to do it in ra-dec, duh)

     if self.InfoFlags[2]=='0':
        L_mean_lat = +1.
     else:
        L_mean_lat = np.mean(self.angular_momentum_helio[1])

     l, b = pole_track_helio.galactic.l, pole_track_helio.galactic.b
     pole_track_helio = ac.SkyCoord(l=l, b=b, frame='galactic')
     #Flip pole track to match Lsign only if it's negative, which can only happen if L exists, if not, it is set to >0 by default
     if L_mean_lat<0 and np.mean(pole_track_helio.galactic.b)>0 :
     	m = b>0.*u.deg  
     	l[m] = l[m] + 180.*u.deg
     	b[m] = -b[m]
     	pole_track_helio = ac.SkyCoord(l=l, b=b, frame='galactic')

     if L_mean_lat>0 and np.mean(pole_track_helio.galactic.b)<0 :  
     	m = b<0.*u.deg
     	l[m] = l[m] + 180.*u.deg
     	b[m] = np.abs(b[m])
     	pole_track_helio = ac.SkyCoord(l=l, b=b, frame='galactic')
        

     #Compute galactocentric pole now. Poles transform as pole(r1,r2)_gc = pole(r1,r2)_helio + (r1-r2)x(rsun_wrt_gc)
     #I can do that, or as done here, trasnform first to GSR and then compute the pole as before
     ep1_gc = ep1.transform_to(ac.Galactocentric())
     ep2_gc = ep2.transform_to(ac.Galactocentric())

     #Will return galactocentric pole as well
     pole_track_gsr = gc.pole_from_endpoints(ep1_gc,ep2_gc).spherical
     lon, lat = pole_track_gsr.lon, pole_track_gsr.lat
#     m = lat<0.*u.deg
#     lon[m] = lon[m] + 180.*u.deg
#     lat[m] = np.abs(lat[m])
     pole_track_gsr = ac.SkyCoord(lon=lon, lat=lat, frame=ac.Galactocentric(), representation_type='spherical')

     return pole_track_helio, pole_track_gsr
   
  def create_sky_polygon_footprint_from_track(self, width=1.*u.deg, phi2_offset=0.*u.deg):

    ''' 
      Create the Polygon Footprint from the celestial track. The polygon is created by shifting the track in phi2 by a given width. 
      Default width for now is 1 deg

    '''  

    #Convert to stream's coordinate frame
    tr = self.track.transform_to(self.stream_frame)

    #Create poly by shifting the track N/S in phi2 by a given angular width
    sort = np.argsort(tr.phi1)
    tr_N = ac.SkyCoord(phi1 = tr.phi1[sort], phi2 = tr.phi2[sort] + width + phi2_offset, frame=self.stream_frame)
    tr_S = ac.SkyCoord(phi1 = tr.phi1[sort], phi2 = tr.phi2[sort] - width + phi2_offset, frame=self.stream_frame)

    #Set poly
    poly_sc = ac.SkyCoord(phi1 = np.append(tr_N.phi1,tr_S.phi1[::-1]) , phi2 = np.append(tr_N.phi2,tr_S.phi2[::-1]), unit=u.deg, frame=self.stream_frame)

    return poly_sc

  def get_mask_in_poly_footprint(self,coo):

     ''' Return a mask for  points in input SkyCoords object that are inside polygon footprint. 

         Parameters
	 ==========

	 coo : astropy.coordinates.SkyCoord object

         Returns
	 =======

	 mask : boolean mask array, same number of elements as coo 
     '''

     return get_mask_in_poly_footprint(poly_sc=self.poly_sc, coo=coo, stream_frame=self.stream_frame)


  def resample_stream_track(self, dphi1=0.02*u.deg): 

     ''' TODO '''


  #-------------method to plot whole MW streams compilation object at once------------------------------------
  def plot_stream_compilation(self,ax,Rstat='mean',Rrange=[0.,9e9],plot_stream_type='all',plot_names=True,
                              use_shortnames=False,plot_colorbar=False,
                              scat_kwargs=None,text_kwargs=None,sym_kwargs=None,cb_kwargs=None,cootype='gal',verbose=False,
                              exclude_streams=[],include_only=[]): 

   '''plot_stream_compilation(self,ax,Rstat='mean',Rrange=[0.,9e9],plot_stream_type='all',plot_names=True,
                              use_shortnames=False,plot_colorbar=False,
                              scat_kwargs=None,text_kwargs=None,sym_kwargs=None,cb_kwargs=None,cootype='gal',verbose=False,
                              exclude_streams=[],include_only=[])  '''


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
