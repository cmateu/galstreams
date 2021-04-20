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


#---------------------------------
#Footprint class definition
class Footprint:
    def __init__(self,lon,lat,name,Dist=None,vrad=None,pmlon=None,pmlat=None,cootype='gal',degree=True,is_pml_star=True,
                 xyz_sun=[8.5,0.,0.],vel_sun=[10.3,232.6,5.9],vxyz_gal=(),Rhel=None):
        
        self.deg=degree 
        self._f=np.pi/180.  
        self.name=name
        self.sname=name[:8]      

        #Bovy's library assumes Sun's position is positive. Flip X-axis is done later
        sign=1.
        #Save Sun's position
        self.xsun, self.ysun, self.zsun= xyz_sun
        self.vxsun, self.vysun, self.vzsun= vel_sun
        self.xsun, self.vxsun = sign*self.xsun, sign*self.vxsun
        
        #Rhel is kept for backwards compatibility
        if Rhel is not None and 'GC' in cootype: 
          sys.exit('ERROR: Rhel not compatible with cootype=GC (distance provided should be galactocentric)')
        elif Rhel is not None: Dist,self.Rhel=Rhel,Rhel

        if 'GC' in cootype:
          if Dist is not None: self.Rgal=Dist
          else: sys.exit('ERROR: Distance is mandatory to instantiate a Footprint in Galactocentric coords (cootype=GC)')

          self.phi,self.theta=lon,lat
          if vxyz_gal: self.vx,self.vy,self.vz=vxyz_gal
          #Compute and set all heliocentric attributes
          self.compute_heliocentric_coords(degree=degree)
        else:
          if 'gal' in cootype:
           self.l,self.b=lon,lat
           mm=bovyc.lb_to_radec(self.l,self.b,degree=self.deg)  
           if self.l.size>1: self.ra,self.dec=mm.T
           else: self.ra,self.dec=mm       
          else:
           self.ra,self.dec=lon,lat 
           mm=bovyc.radec_to_lb(self.ra,self.dec,degree=self.deg)
           if self.ra.size>1: self.l,self.b=mm.T 
           else: self.l,self.b=mm
              
          if Dist is not None: self.Rhel=Dist
          if vrad is not None: self.vrad=vrad
            
          if 'gal' in cootype:
            if pmlon is not None: 
                if is_pml_star: self.pmlstar,self.pml=pmlon,pmlon/np.cos(_f*self.b)
                else: self.pmlstar,self.pml=pmlon*np.cos(self._f*self.b),pmlon
            if pmlat is not None: self.pmb=pmlat
          else:        
            if pmlon is not None: 
                if is_pml_star: self.pmrastar,self.pmra=pmlon,pmlon/np.cos(_f*self.dec)
                else: self.pmrastar,self.pmra=pmlon*np.cos(self._f*self.dec),pmlon
            if pmlat is not None: self.pmdec=pmlat

          #Set galactocentric attributes
          if hasattr(self,'Rhel') : self.compute_galactocentric_coords(degree=degree)

        #Set center attributes
        self.compute_sky_center()

        #Compute and Set end-point attributes
        try:
           self.compute_midplane_endpoints_1(verbose=False) 
           self.mp=1
        except:
           self.compute_midplane_endpoints_2(verbose=False) 
           self.mp=2

        #Set decent astropy Skycoord object as attribute
        self.sc = astropy.coordinates.SkyCoord(ra=self.ra*u.deg,dec=self.dec*u.deg)

        #Set great-circle gala-reference-frame for each stream based on its mid-plane end-points
        self.gcfr = gc.GreatCircleICRSFrame.from_endpoints(self.end_o, self.end_f)

        #Check that pole is not nan, use method_1 to correct it if it is, and recompute gcfr
        if not (self.gcfr.pole.ra>=0):
           self.compute_midplane_endpoints_2(verbose=False,tol=0.1) 
           self.mp=2
           self.gcfr = gc.GreatCircleICRSFrame.from_endpoints(self.end_o, self.end_f)
        #Flip if pole's dec is negative
        if self.gcfr.pole.dec<0:
          self.gcfr=gc.GreatCircleICRSFrame.from_endpoints(self.end_f,self.end_o) 
 
        #Provide phi1 and phi2 as "normal" Footprint attributes
        self.phi1 = self.sc.transform_to(self.gcfr).phi1  
        self.phi2 = self.sc.transform_to(self.gcfr).phi2  
       
        #goes here

      
    def compute_sky_center(self):           
        
        #Set defaults        
        #Need to get cartesian coords to do vector-average (this works whether or not Rhel exists)
        mm=bovyc.lbd_to_XYZ(self.l,self.b,np.ones_like(self.l),degree=True)
        if self.l.size>1: _xx,_yy,_zz=mm.T
        else: _xx,_yy,_zz=mm    
        _xc,_yc,_zc=_xx.sum(),_yy.sum(),_zz.sum()
        self.cl,self.cb=bovyc.XYZ_to_lbd(_xc,_yc,_zc,degree=True)[:2]
        self.cra,self.cdec=bovyc.lb_to_radec(self.cl,self.cb,degree=True)

        if hasattr(self,'phi') and hasattr(self,'theta'):
          _xgc,_ygc,_zgc=self.x.sum(),self.y.sum(),self.z.sum()
          self.cphi,self.ctheta=bovyc.XYZ_to_lbd(_xgc,_ygc,_zgc,degree=True)[:2]


    def mask_footprint(self,mask):
        
        #Apply mask to all array attributes of object
        for myattr in self.__dict__.keys():
        #If attribute is callable, proceed (only applicable to ndim>=1 arrays)
          if not callable(getattr(self,myattr)) and np.ndim(getattr(self,myattr))>=1:
            setattr(self,myattr,getattr(self,myattr)[mask])
            #print 'Attribute ',myattr, len(getattr(selfcopy,myattr))


#---------MW Streams class--------------------------------------------------------------------------------
        
class MWStreams(dict):
    
  def __init__(self,verbose=True):

    #A MWStreams object is a dictionary in which each entry is a Footprint object, indexed by each stream's name.
    #There's also a mandatory summary entry, which is a Pandas DataFrame with summary attributes for the full library

    #Initialize empty dictionary
    dict.__init__(self)

    #Read in the master log
    lib_master_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),'lib','master_log.txt')

    tdir = os.path.dirname(os.path.realpath(__file__))
    lmaster = astropy.table.Table.read(tdir+"/lib/master_log.txt", format='ascii.commented_header').to_pandas()
    lmaster["On"] = lmaster["On"].astype('bool') #this attribute controls whether a given track is included or not

    #SkyCoords objects will be created for each of these dicts after the full library has been created
    attributes = ['ra','dec','distance','pm_ra_cosdec','pm_dec','radial_velocity']
    units = [u.deg, u.deg, u.kpc, u.mas/u.yr, u.mas/u.yr, u.km/u.s]
    end_o_dic = {k: np.array([])*uu  for k,uu in zip(attributes,units) }
    end_f_dic = {k: np.array([])*uu  for k,uu in zip(attributes,units) }
    mid_point_dic = {k: np.array([])*uu  for k,uu in zip(attributes,units) }
    mid_pole_dic  = {k: np.array([])*uu  for k,uu in zip(attributes,units) }
    info_flags = []

    print("Initializing galstreams library from master_log... ")
    for ii in np.arange(lmaster.Name.size):

       print(lmaster.Name[ii])

       #Create the names of the files containing the knots and summary attributes to initialize each stream
       summary_file = "{tdir}/track.{imp}.{stname}.{tref}.summary.ecsv".format(tdir=tdir+'/tracks', imp=lmaster.Imp[ii], 
          								      stname=lmaster.Name[ii],
          								      tref=lmaster.TrackRefs[ii])
       track_file   = "{tdir}/track.{imp}.{stname}.{tref}.ecsv".format(tdir=tdir+'/tracks', imp=lmaster.Imp[ii], 
          							      stname=lmaster.Name[ii],
          							      tref=lmaster.TrackRefs[ii])

       if verbose: 
          print(f"Initializing Track6D {TrackName[ii]}...")

       #Do the magic. The track is read and all attributes stored in the summary for all registered stream tracks. 
       #But only the ones "On" are "realized"
       track = Track6D(track_name=lmaster.TrackName[ii], stream_name=lmaster.Name[ii], track_file=track_file, summary_file=summary_file)

       if lmaster.On[ii]: self[lmaster.TrackName[ii]] = track
       else: 
          if verbose: print(f"Skipping {TrackName[ii]}...")

       #Store summary attributes
       for k in attributes: end_o_dic[k] = np.append(end_o_dic[k], getattr(track.end_points, k)[0] )
       for k in attributes: end_f_dic[k] = np.append(end_f_dic[k], getattr(track.end_points, k)[1] )
       for k in attributes: mid_point_dic[k] = np.append(mid_point_dic[k], getattr(track.mid_point, k) )
       for k in attributes[:2]: mid_pole_dic[k]  = np.append(mid_pole_dic[k] , getattr(track.mid_pole, k) )

       info_flags.append(track.info_flags)

    #Store master table as an attribute
    self.summary = lmaster

    #Add skycoord summary attributes to the library and selected cols to the summary table
    self.end_o = ac.SkyCoord(**end_o_dic)
    self.end_f = ac.SkyCoord(**end_f_dic)
    self.mid_point = ac.SkyCoord(**mid_point_dic)
    self.mid_pole  = ac.SkyCoord(ra=mid_pole_dic["ra"], dec=mid_pole_dic["dec"], frame='icrs')

    self.summary["ra_o"] = end_o_dic["ra"]
    self.summary["dec_o"] = end_o_dic["dec"]
    self.summary["ra_f"] = end_f_dic["ra"]
    self.summary["dec_f"] = end_f_dic["dec"]
    self.summary["ra_mid"] = mid_point_dic["ra"]
    self.summary["dec_mid"] = mid_point_dic["dec"]
    self.summary["distance"] = mid_point_dic["distance"]
    self.summary["ra_pole"] = mid_pole_dic["ra"]
    self.summary["dec_pole"] = mid_pole_dic["dec"]
    self.summary["InfoFlags"] = np.array(info_flags)


class Track6D:

  def __init__(self, track_name, track_file, summary_file, stream_name=None, stream_shortname=None, verbose=True):

      ''' Track6D: A Stellar Stream's Track realization in 6D. The object has the following attributes: 
 	
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

	end_points: 2-element astropy.coordinates.SkyCoord Object
	  
        mid_point: astropy.coordinates.SkyCoord Object

        mid_pole: astropy.coordinates.SkyCoord Object

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
      else: self.stream_shortname = self.stream_name[:3]

      #Read-in knots and initialize track
      t = astropy.table.QTable.read(track_file)

      #Store the track in attribute
      self.track = ac.SkyCoord(**t)

      #Now read in the stream's summary file 
      sfile=astropy.table.QTable.read(summary_file)

      #Streams InfoFlags: four-bit flag
      # bit 0: 0 = great circle by construction
      # bit 1: 0 = no distance track available (only mean or central value reported)
      # bit 2: 0 = no proper motion track available (only mean or central value reported)
      # bit 3: 0 = no radial velocity track available (only mean or central value reported)
      self.info_flags = str(sfile["InfoFlags"][0]) # All-in-one flag

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
      x["distance"] = 1.*u.kpc   #it shouldn't matter, but if it's zero it does
      self.mid_pole = ac.SkyCoord(**x)

      #Set up stream's coordinate frame
      self.stream_frame = gc.GreatCircleICRSFrame(pole=self.mid_pole, ra0=self.mid_point.icrs.ra)

      self.poly_sc = self.create_sky_polygon_footprint_from_track(width=1*u.deg)

  def compute_pole_track(self):

     ''' TODO '''
   
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

     ''' Test whether points in input SkyCoords object are inside polygon footprint. 

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


  def load_user_defined_centers_and_shortnames(self):

     #Read library log-file that will be used to overwrite center coords with user-defined values
     lib_log_filen=os.path.join(os.path.dirname(os.path.realpath(__file__)),'lib','lib_centers.log')
     names,shortnames=scipy.genfromtxt(lib_log_filen,usecols=(0,1),unpack=True,dtype=str)
     _ra,_dec,_ll,_bb,_phi,_theta=scipy.genfromtxt(lib_log_filen,usecols=(2,3,4,5,6,7),
                                                    unpack=True,dtype=np.float,missing_values='',filling_values=np.nan)
 
     for ii in range(names.size):
       try:
         if not np.isnan( _ra[ii])  and not np.isnan(  _dec[ii]): self[names[ii]].cra, self[names[ii]].cdec =  _ra[ii],   _dec[ii] 
         if not np.isnan( _ll[ii])  and not np.isnan(   _bb[ii]):  self[names[ii]].cl, self[names[ii]].cb   =  _ll[ii],    _bb[ii] 
         if  hasattr(self[names[ii]],'cphi'):
          if not np.isnan(_phi[ii])  and not np.isnan(_theta[ii]): self[names[ii]].cphi,self[names[ii]].ctheta = _phi[ii], _theta[ii] 
         #Setting short names
         self[names[ii]].sname=shortnames[ii]        
       except KeyError: print('WARNING: Name %s used lib_centers.log not found in lib* definition files' % (names[ii]))


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
