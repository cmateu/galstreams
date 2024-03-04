import numpy as np
import scipy
import matplotlib as mpl
import pylab as plt
import os, os.path
import sys
import astropy.table
import astropy.coordinates as ac
import astropy.units as u
import gala
import gala.coordinates as gc
import gala.dynamics as gd
import pandas as pd

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

def create_sky_polygon_footprint_from_track(SkyCoordTrack, frame, width=1.*u.deg, phi2_offset=0.*u.deg):

  ''' 
    Create the Polygon Footprint from the celestial track. The polygon is created by shifting the track in phi2 by a given width. 
    
    Inputs:
    =======

    SkyCoordTrack: track SkyCoord object from a MWStreams library stream (mws[st].track)

    frame: None. Astropy coordinate frame to set up the polygon by offsetting the track by a given width. 
           The default is to use the Track6D's own stream frame Track6D.stream_frame

    Parameters
    ==========

    phi2_offset: astropy.Quantity object
     The offset in phi2 that will be applied to the track to create the polygon footprint (default 0)

    width: astropy.Quantity object
     The total width of the polygon footprint to be created around track+phi2_offset

  '''

  track = SkyCoordTrack
  #if frame is None: 
  # frame = Track6D.stream_frame

  #Convert to stream's coordinate frame
  tr = track.transform_to(frame)

  #Create poly by shifting the track N/S in phi2 by a given angular width
  sort = np.argsort(tr.phi1)
  tr_N = ac.SkyCoord(phi1 = tr.phi1[sort], phi2 = tr.phi2[sort] + width/2. + phi2_offset, frame=frame)
  tr_S = ac.SkyCoord(phi1 = tr.phi1[sort], phi2 = tr.phi2[sort] - width/2. + phi2_offset, frame=frame)

  #Set poly
  # Concatenate N track, S-flipped track and add first point at the end to close the polygon (needed for ADQL)
  poly_sc = ac.SkyCoord(phi1 = np.concatenate((tr_N.phi1,tr_S.phi1[::-1],tr_N.phi1[:1])),
                        phi2 = np.concatenate((tr_N.phi2,tr_S.phi2[::-1],tr_N.phi2[:1])),
                        unit=u.deg, frame=frame)

  return poly_sc


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
    
  def __init__(self, verbose=False, implement_Off=False, print_topcat_friendly_files=False):

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
    #separate info_flags (easier to filter)
    has_empirical_track = np.array([],dtype=np.int32)
    has_D = np.array([],dtype=np.int32)
    has_pm = np.array([],dtype=np.int32)
    has_vrad = np.array([],dtype=np.int32)
    discovery_refs = []
    lengths = np.array([])#*u.deg

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
       has_empirical_track = np.append(has_empirical_track, np.int32(track.InfoFlags[0]))
       has_D               = np.append(has_D    , np.int32(track.InfoFlags[1])) 
       has_pm              = np.append(has_pm   , np.int32(track.InfoFlags[2]))  
       has_vrad            = np.append(has_vrad , np.int32(track.InfoFlags[3]))  
       lengths = np.append(lengths, track.length.deg)
       discovery_refs = np.append(discovery_refs, lmaster_discovery.loc[lmaster.Name[ii],'DiscoveryRefs'] )


    #Add skycoord summary attributes to the library and selected cols to the summary table
    self.end_o = ac.SkyCoord(**end_o_dic)
    self.end_f = ac.SkyCoord(**end_f_dic)
    self.mid_point = ac.SkyCoord(**mid_point_dic)
    self.mid_pole  = ac.SkyCoord(ra=mid_pole_dic["ra"], dec=mid_pole_dic["dec"], frame='icrs')

    #Store master table as an attribute (inherits structure of lmaster dataframe)
    self.summary = lmaster.copy()

    #Stream Length
    self.summary["length"] = np.array(lengths)
    #End points
    self.summary["ra_o"] = end_o_dic["ra"].deg
    self.summary["dec_o"] = end_o_dic["dec"].deg
    self.summary["distance_o"] = end_o_dic["distance"].value
    self.summary["ra_f"] = end_f_dic["ra"].deg
    self.summary["dec_f"] = end_f_dic["dec"].deg
    self.summary["distance_f"] = end_f_dic["distance"].value
    #Mid point
    self.summary["ra_mid"] = mid_point_dic["ra"].deg
    self.summary["dec_mid"] = mid_point_dic["dec"].deg
    self.summary["distance_mid"] = mid_point_dic["distance"].value
    #Pole
    self.summary["ra_pole"] = mid_pole_dic["ra"].deg
    self.summary["dec_pole"] = mid_pole_dic["dec"].deg
    #Info (InfoFlags and has_* columns is the same, but to have it on separate columns is more practical for filtering)
    self.summary["InfoFlags"] = np.array(info_flags)
    self.summary["has_empirical_track"] = has_empirical_track
    self.summary["has_D"]    = has_D              
    self.summary["has_pm"]   = has_pm             
    self.summary["has_vrad"] = has_vrad           
    self.summary["DiscoveryRefs"] = discovery_refs

    #Index by TrackName
    self.summary.index=self.summary.TrackName

    #Create a numeric ID for each track
    self.summary["ID"] = ''
    for ii in self.summary.index:
       if self.summary.loc[ii,'On']:
        self.summary.loc[ii,'ID'] = self[ii].ID
  
    #If chosen by the user, when the library is instantiated, save in default location TOPCAT-friendly csv files with 
    # the library's tracks, end-points, mid-points and summary table 
    if print_topcat_friendly_files:
     self.print_topcat_friendly_compilation(output_root=f'{tdir}/tracks/galtreams.unique_streams')


  def all_unique_stream_names(self):
       '''
         Returns all unique instances of the StreamNames in the library (a stream can have multiple tracks)

         Returns
         =======
         
         array 
       '''
       return np.unique(self.summary.Name[self.summary.On])

  def all_track_names(self, On_only=False):
       '''
         Returns TrackNames available in the library (when On_only=False, equivalent to MWStreams.summary['TrackName']) 
  
         Parameters:
         ===========
     
         On_only: True/False
                  If True it returns only the names for the active tracks

         Returns
         =======
         
         array  
       '''
 
       if On_only: return self.keys()  
       else:
             return np.array(self.summary.index)       

  def get_track_names_for_stream(self, StreamName, On_only=False):
    '''
        Find all the TrackNames for which the StreamName matches the input string (all or part of it)

        Parameters
        ==========

        StreamName : str 
                     Name of the stream for which to search TrackNames (or part of it)
        On : book
             If True, returns only the active track name(s)	
     
        Returns 
        =======

        array : contains all the TrackNames for which the stream's name matches the input string  

    '''

    all_track_names = []

    if On_only: track_names = self.keys()
    else: track_names = self.summary.index

    for tn in track_names:
      if StreamName in self.summary.loc[tn,'Name']: 
         all_track_names.append(tn)

    return all_track_names

  def print_topcat_friendly_compilation(self, output_root='galtreams.unique_streams'):

        modes = ['track','end_point','mid_point']
	
        for mode in modes:
            
            full = dict(ra=np.array([]), dec=np.array([]), distance=np.array([]), pm_ra_cosdec=np.array([]), pm_dec=np.array([]), 
                    radial_velocity=np.array([]), ID=np.array([]), StreamName=np.array([]), TrackName=np.array([]))
            
            for st in np.sort(list(self.keys())):
             if self.summary.loc[st,"On"]:
            
                if 'track' in mode: x = self[st].track
                elif 'end' in mode: x = self[st].end_points
                elif 'mid' in mode: x = self[st].mid_point
                N=np.size(x.ra)
            
                full['ra']  = np.append(full['ra'], x.ra.deg)
                full['dec'] = np.append(full['dec'], x.dec.deg)
                full['distance'] = np.append(full['distance'], x.distance.value)
                full['pm_ra_cosdec'] = np.append(full['pm_ra_cosdec'], x.pm_ra_cosdec.value)
                full['pm_dec'] = np.append(full['pm_dec'], x.pm_dec.value)
                full['radial_velocity'] = np.append(full['radial_velocity'], x.radial_velocity.value)
                full['ID'] = np.append(full['ID'], self.summary.loc[st,"ID"] + np.zeros(N, dtype=int))
                full['TrackName'] = np.append(full['TrackName'], [self[st].track_name,]*N)
                full['StreamName']= np.append(full['StreamName'], [self[st].stream_name,]*N)
            
            full_pd = pd.DataFrame.from_dict(full)
            print(f"Creating TOPCAT-friendly library files:\n {output_root}.tracks/end_points/mid_points.csv")
            full_pd.to_csv(f'{output_root}.{mode}s.csv')    
        #Print summary table
        print(f" {output_root}.summary.csv")
        self.summary.to_csv(f'{output_root}.summary.csv')

  def get_track_names_in_sky_window(self, lon_range, lat_range, frame, On_only=True, wrap_angle=0.*u.deg):
    ''' 
       Get a list of track names for streams in a sky window with limits given by 
       lon_range,lat_range in a given coordinate frame 

       Parameters
       ==========

       lon_range : np.array
                   2-element array containing limits of sky window in "longitude" coordinate (e.g ra, l)                  
 
       lat_range : np.array
                   2-element array containing limits of sky window in "latitude" coordinate (e.g dec, b)                  
       
       frame : AstroPy coordinate frame
               Coordinate frame corresponding to lon/lat_range coordinates provided above
  
    '''
        
    #This is just so I can get the representation_component_names (don't know how to do it 
    #without creating a frame instance, so, there, let's move on
    coo = ac.SkyCoord(lon_range, lat_range,frame=frame, unit=u.deg)
    n = dict((v,k) for k,v in coo.frame.representation_component_names.items())
        
    if np.any(lon_range<0.): 
        wrap_angle = 180.*u.deg
        #print(f"Warning: negative longitudes - setting wrap_angle to 180. assuming this is not a mistake")
        
    llon = ac.Angle(lon_range).wrap_at(wrap_angle).deg
    llat = ac.Angle(lat_range).deg

    track_names = [] 
    for st in self.summary.TrackName:
        if On_only:
            if ~self.summary["On"][st]: continue
        #same for current track  
        lon = getattr(self[st].track.transform_to(frame), n['lon']).wrap_at(wrap_angle).deg
        lat = getattr(self[st].track.transform_to(frame), n['lat']).deg

        mask =  (np.min(llon)<= lon) & (lon <= np.max(llon)) & (np.min(llat)<= lat) & (lat <= np.max(llat))

        if (mask.sum()>0): track_names.append(st)

    return track_names

  def plot_stream_compilation(self, ax=None, frame=ac.ICRS, C_attribute=None, plot_names='ID',
                              plot_colorbar = None, invert_axis=True, show_legend=True,
                              basemap_kwds = dict(projection='moll',lon_0=180., resolution='l'), 
                              mlabels_kwds = dict(meridians=np.arange(0.,360.,30.), color=(0.65,0.65,0.65),linewidth=1., laxmax=90.),
                              plabels_kwds = dict(circles=np.arange(-75,75,15.), color=(0.65,0.65,0.65),linewidth=1.,
                                                  labels=[0,1,1,0], labelstyle='+/-' ),
                              scat_kwds = None,
                              annot_kwds = dict(xytext=(15,15),
                                                textcoords='offset points',
                                                arrowprops=dict(arrowstyle="-",color='k'),
                                                horizontalalignment='center', verticalalignment='center'), 
			      legend_kwds = dict(ncol=8,loc='center', columnspacing=0.5, handletextpad=0.1,
 			                         bbox_to_anchor=(0.5,-0.28), markerscale=3, fontsize='medium'),
                              cb_kwds = None,
                              exclude_streams=[], include_only=[], 
                              return_basemap_m = False,
  			      verbose=False): 
   ''' 

     Plot a Mollweide sky projection map of the current MWStreams library object in the selected coordinate frame.
     Note: requires Basemap library

     Parameters
     ==========

             track : SkyCoord object


         ax=None 

         frame : Astropy astropy.coordinates.baseframe instance
                 Coordinate frame to be used in sky plot

         C_attribute : name of SkyCoord object attribute (in selected reference frame) 
                       to pass plt.scatter as auxiliary column c  
                       e.g. 'distance', 'pm_b' if frame=ac.Galactic

         plot_names : str ['ID','track_name','stream_name','stream_shortname']
         
         plot_colorbar: Bool
                        If C_attribute is passed, plot_colorbar=True by default	

         invert_axis : Bool
                       Invert longitude axis, set to True by default to follow usual plotting convention for l/ra
         
         show_legend: Bool
                      Show legend at the bottom of the plot. Legend attributes can be passed via the legend_kwds dict
         
         basemap_kwds : dict
                        Keywords to instantiate Basemap projection. Default, Molweide projection

         mlabels_kwds: dict  - default=dict(meridians=np.arange(0.,360.,30.), color=(0.65,0.65,0.65),linewidth=1., laxmax=90.)
                       Meridian labelling keyword attributes to be passed to Basemap

         plabels_kwds: dict - default=dict(circles=np.arange(-75,75,15.), color=(0.65,0.65,0.65),linewidth=1.,
                          labels=[0,1,1,0], labelstyle='+/-' )
                       Parallel labelling keyword attributes to be passed to Basemap

         scat_kwds : dict - default scat_kwds=dict(marker='.', s=30, alpha=0.8) [defaults change if C_attribute is passed]
                     Plotting keyword attributes to be passed to plt.scatter

         annot_kwds : dict
                      Text and arrow attributes to be passed to annotate

         legend_kwds : dict
                       Legend attributes to be passed to plt.legend

         cb_kwds : dict - default = dict(label=C_attribute,  shrink=0.5)
                   Colorbar attributes to be passed to plt.colorbar 
         
         exclude_streams: list of stream TrackNames
                          TrackNames for streams *not* to be included in the plot        

         include_only: list of stream TrackNames
                       Only the TrackNames provided in this list will be plotted

         return_basemap_m : False
                            Return Basemap projection function 
 
         verbose: False
                  Not doing anything right now if set to True

     Returns
     =======

     	ax

     	ax : Current axes object

   '''

   if ax is None:
     fig = plt.figure(1,figsize=(16.5,11))
     ax = fig.add_subplot(111)
   
   #Follow the usual convention for Galactic and ICRS to invert the l/ra axis
   if invert_axis:
     ax.invert_xaxis()
   
   if 'Basemap' not in sys.modules:
     from mpl_toolkits.basemap import Basemap
     
   m = Basemap(**basemap_kwds)
   
   m.drawmeridians(**mlabels_kwds)
   m.drawparallels(**plabels_kwds)
   m.drawmapboundary()
   
   fr = frame

   if len(include_only)>0 : keys_to_plot = include_only
   else: keys_to_plot = self.keys()
   
   msg_flag = 0

   for st in keys_to_plot:

    if st in exclude_streams: continue

    if self.summary.loc[st,"On"]:

     #Get representation names for selected frame and flip dict around 
     l = self[st].track.transform_to(fr).representation_component_names
     n = dict((v,k) for k,v in l.items())
 
     if plot_names is None: 
       label, alabel = None, None
     elif 'ID' not in plot_names:     
      label = "{Name}".format(Name = getattr(self[st],plot_names) )
      alabel = label
     else:
         label="{ID:.0f}={Name}".format(ID=self[st].ID,Name=self[st].track_name)
         alabel="{ID:.0f}".format(ID=self[st].ID)
      

     #Transform the current stream's track to selected frame
     coo = self[st].track.transform_to(fr)
     x,y = m( getattr(coo,n['lon']).value , getattr(coo, n['lat']).value )

     if scat_kwds is None: 
      if C_attribute is None:
       scat_kwds=dict(marker='.', s=30, alpha=0.8)
      elif C_attribute == 'distance':
        scat_kwds=dict(marker='.', s=30, alpha=0.8, vmin=0., vmax=100.) #reasonable limits for distance plot
      else:
        scat_kwds=dict(marker='.', s=30, alpha=0.8, vmin=-10., vmax=10.) #reasonable limits to plot pms

     #Extra attribute to plot
     if C_attribute is not None: 
       try: c = getattr(coo, C_attribute).value
       except AttributeError:
          c = None
          if msg_flag==0:
            print('WARNING: Invalid attribute selected. If not a spelling error, you are probably trying to plot an attribute in a different coord frame as the one selected. This is currently not supported. Plotting without C_attribute aux column for now...')
          msg_flag = 1
     else: c = None

     im = ax.scatter(x,y, c=c,  **scat_kwds, label=label)


     #Using end_point to place labels
     coo = self[st].mid_point.transform_to(fr)
     xl,yl = m( getattr(coo,n['lon']).value , getattr(coo, n['lat']).value )
     xy_stream = xl, yl  
     ax.annotate(alabel, xy=xy_stream,  xycoords='data', **annot_kwds)

   ax.grid(ls=':')

   if show_legend:
    ax.legend(**legend_kwds)

   if cb_kwds is None and C_attribute is not None: 
     cb_kwds = dict(label=C_attribute,  shrink=0.5)

   if C_attribute is not None and plot_colorbar is None: plot_colorbar=True 
   if plot_colorbar:
     plt.colorbar(im, ax=ax, **cb_kwds)
 
   if return_basemap_m:
      return ax, m
   else:
      return ax

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

	stream_frame: astropy coordinate frame  

	track : astropy.coordinates.SkyCoord Object
	  Contains the track 6D info. By default initialized in icrs frame

        length: astropy.Quantity Object contains the angular length measured along the track

        InfoFlags: string - 4-bits indicate available (or assumed) data. 
  	    bit 0: 0 = great circle by construction
  	    bit 1: 0 = no distance track available (only mean or central value reported)
  	    bit 2: 0 = no proper motion data available (only mean or central value reported)
  	    bit 3: 0 = no radial velocity data available (only mean or central value reported)
 
	end_points: 2-element astropy.coordinates.SkyCoord Object with end point coordinates
	  
        mid_point: astropy.coordinates.SkyCoord Object with stream's mid-point coordinates (phi1=0)

        mid_pole: astropy.coordinates.SkyCoord Object heliocentric pole at mid_point
 
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
#<<<<<<< patch-1
#      self.stream_frame = gc.GreatCircleICRSFrame.from_pole_ra0(pole=self.mid_pole, ra0=self.mid_point.icrs.ra)
#=======
      #for now
      if np.float64(gala.__version__[:3])<=1.4:
         self.stream_frame = gc.GreatCircleICRSFrame(pole=self.mid_pole, ra0=self.mid_point.icrs.ra)  
      else:  
         self.stream_frame = gc.GreatCircleICRSFrame.from_pole_ra0(
            pole=self.mid_pole,
            ra0=self.mid_point.icrs.ra,
            origin_disambiguate=self.mid_point.icrs
         )
#>>>>>>> master

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

     ''' Compute pole at each point in the track. This is obtained by computing, at each point, the normal or cross product betwee
         said point and the contiguous point in the track '''

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

    poly_sc = create_sky_polygon_footprint_from_track(self.track, frame=self.stream_frame, width=width, phi2_offset=phi2_offset)

    return poly_sc

  def get_adql_query_from_polygon(self):

     return get_adql_query_from_polygon(self.poly_sc)

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

     ''' In construction... '''

   

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
