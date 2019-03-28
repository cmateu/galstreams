import numpy as np
import scipy
import pylab as plt
import bovy_coords as bovyc
import gcutils 
import os, os.path
import sys
import astropy.coordinates
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

    def compute_galactocentric_coords(self,verbose=False,degree=True):
      #Convert to heliocentric cartesian. Bovy's library assumes Sun's position is positive
      #tuple output, no .T needed
      if verbose: print('Converting Heliocentric Galactic Spherical to Heliocentric Cartesian coords...')
      m=bovyc.lbd_to_XYZ(self.l,self.b,self.Rhel,degree=degree)
      self.xhel,self.yhel,self.zhel=m.T
      if hasattr(self,'vrad') and hasattr(self,'pmb') :
        m=bovyc.vrpmllpmbb_to_vxvyvz(self.vrad,self.pmlstar,self.pmb,self.l,self.b,self.Rhel,XYZ=False,degree=degree)
        self.vxhel,self.vyhel,self.vzhel=m.T

      #Convert Heliocentric Cartesian to Galactocentric Cartesian
      if verbose: print('Converting Heliocentric Cartesian to Galactocentric Cartesian coords...')
      m=bovyc.XYZ_to_galcenrect(self.xhel,self.yhel,self.zhel,Xsun=self.xsun,Ysun=self.ysun,Zsun=self.zsun)
      self.x,self.y,self.z=m
      if hasattr(self,'vrad') and hasattr(self,'mub') :
       m=bovyc.vxvyvz_to_galcenrect(self.vxhel,self.vyhel,self.vzhel,vsun=[self.vxsun, self.vysun, self.vzsun])
       self.vx,self.vy,self.vz=m
  
      #Compute Galactocentric Spherical
      m=bovyc.XYZ_to_lbd(self.x,self.y,self.z,degree=degree)
      self.phi,self.theta,self.Rgal=m.T

      #Bovy's ref system's x-axis points towards the Sun. Mine points away from the Sun, i.e. x-axis is the same as for lbd
      #Will use my ref system for consistency with pole-ref system used in PyMGC3
      if degree: f=180./np.pi
      else: f=1.
      self.x=-self.x
      self.phi= (2*np.pi*f - self.phi) - np.pi*f
      self.phi[self.phi<0]=self.phi[self.phi<0]+2*np.pi*f

    def compute_heliocentric_coords(self,verbose=False,degree=True):

       #Set Galactocentric Cartesian coords
       xg,yg,zg=bovyc.lbd_to_XYZ(self.phi,self.theta,self.Rgal,degree=degree).T
       self.x,self.y,self.z=-xg,yg,zg
       #Convert to Heliocentric
       self.xhel,self.yhel,self.zhel=bovyc.galcenrect_to_XYZ(self.x,self.y,self.z,
                                                             Xsun=self.xsun,Ysun=self.ysun,Zsun=self.zsun)
       #Save Heliocentric galactic and equatorial
       self.l,self.b,self.Rhel=bovyc.XYZ_to_lbd(self.xhel,self.yhel,self.zhel,degree=degree).T
       self.ra,self.dec=bovyc.lb_to_radec(self.l,self.b,degree=degree).T
 
       #Set kinematic attributes, if vels available   
       if hasattr(self,'vx') and hasattr(self,'vy') and hasattr(self,'vz'):
          #Cartesian Heliocentric
          m=bovyc.galcenrect_to_vxvyvz(self.vx,self.vy,self.vz,vsun=[self.vxsun, self.vysun, self.vzsun])
          self.vxhel,self.vyhel,self.vzhel=m
          #Spherical Heliocentric Galactic
          m=bovyc.vxvyvz_to_vrpmllpmbb(self.vxhel,self.vyhel,self.vzhel,self.l,self.b,self.Rhel,XYZ=False,degree=degree)
          self.vrad,self.pmlstar,self.pmb=m.T
          self.pml=self.pmlstar/np.cos(self.b*self._f)  
          #Spherical Heliocentric Equatorial
          self.pmrastar,self.pmdec=pmllpmbb_to_pmrapmdec(self.pml,self.pmb,self.l,self.b,degree=degree,epoch=2000.0) 
          self.pmra=self.pmrastar/np.cos(self.dec*self._f)  

    def compute_midplane_endpoints_2(self,verbose=False, tol=0.5):

          if verbose: print("Computing end points for %s..." % (self.name))

          ramin, ramax = np.min(self.ra), np.max(self.ra)
          rari_mask = (self.ra<=ramin+tol)
          rale_mask = (self.ra>=ramax-tol)
        
          ra_o,dec_o = get_avg_vec(self.ra[rari_mask],self.dec[rari_mask])
          ra_f,dec_f = get_avg_vec(self.ra[rale_mask],self.dec[rale_mask])            

          #Set the skycoord end-points as attributes
          self.end_o = astropy.coordinates.SkyCoord(ra_o, dec_o,frame='icrs', unit=u.deg)
          self.end_f = astropy.coordinates.SkyCoord(ra_f, dec_f,frame='icrs', unit=u.deg)

    def compute_midplane_endpoints_1(self,verbose=False):

          if verbose: print("Computing end points for %s..." % (self.name))

          ra1,dec1 = self.ra[np.argmin(self.dec)], self.dec[np.argmin(self.dec)]
          ra2,dec2 = self.ra[np.argmin(self.ra)],  self.dec[np.argmin(self.ra)]
          ra3,dec3 = self.ra[np.argmax(self.ra)],  self.dec[np.argmax(self.ra)]
          ra4,dec4 = self.ra[np.argmax(self.dec)], self.dec[np.argmax(self.dec)]

          call = astropy.coordinates.SkyCoord(ra=np.array([ra1,ra2,ra3,ra4])*u.deg,
                                              dec=np.array([dec1,dec2,dec3,dec4])*u.deg)

          triup = np.triu(call[:, None].separation(call[None]).degree)
          idx1, idx2 = np.where(np.isclose(triup, triup[triup != 0].min()))

          end1a = call[idx1[0]]
          end1b = call[idx2[0]]
          end2a = call[idx1[1]]
          end2b = call[idx2[1]]

          #Set the skycoord end-points as attributes
          self.end_o = astropy.coordinates.SkyCoord(gc.greatcircle.sph_midpoint(end1a, end1b),frame='icrs')
          self.end_f = astropy.coordinates.SkyCoord(gc.greatcircle.sph_midpoint(end2a, end2b),frame='icrs')


#---------MW Streams class--------------------------------------------------------------------------------
        
class MWStreams(dict):
    
  #by pair of end-points------------------------------------------------------------------------------
  def init_by_end_points(self, gcstep=0.01,verbose=False):
    
      lib_end_points_filen=os.path.join(os.path.dirname(os.path.realpath(__file__)),'lib','lib_by_pair.dat')

      lono,lato,lonf,latf,ro,rf,width=scipy.genfromtxt(lib_end_points_filen,usecols=(3-1,4-1,5-1,6-1,8-1,9-1,10-1),unpack=True)
      name,sname,cootype=scipy.genfromtxt(lib_end_points_filen,usecols=(1-1,2-1,7-1),unpack=True,dtype=str)
  
      for i in range(len(lono)):
        #Get great-circle lons,lats given end-points 
        if verbose: print('Reading pair-list for %s' % (name[i]))
        azs,lats,azcenter,latcenter=gcutils.get_gc_for_pair(lono[i],lato[i],lonf[i],latf[i],degree=True,step=gcstep,dlat=width[i]) 
        
        #Do linear interpolation for the distance, not much better to do
        if ro[i]==rf[i]: D=ro[i]*np.ones_like(azs)
        elif rf[i]>ro[i]>0.: D=scipy.interp(azs,[lono[i],lonf[i]],[ro[i],rf[i]])    
        else: D=-1*np.ones_like(azs)
        
        #Create footprint appropriately depending on coord type, this will define the l,b,ra,dec attributes
        footprint=Footprint(azs,lats,name[i],Dist=D,degree=True,cootype=cootype[i])
        
        #Store
        self[name[i]]=footprint     
    
  #--------------by pole----------------------------------------------------------
  def init_by_poles(self,gcstep=0.1,verbose=False):  

        lib_poles_filen=os.path.join(os.path.dirname(os.path.realpath(__file__)),'lib','lib_by_pole.dat')

        name,pole_coot,c_coot=scipy.genfromtxt(lib_poles_filen,usecols=(0,3,6),unpack=True,dtype=str)
        pole_dat=scipy.genfromtxt(lib_poles_filen,usecols=(1,2,4,5,7,8,9,10),filling_values=-999.)
        pole_lons,pole_lats,clons,clats,dlons,dlats,ros,rfs=pole_dat.T

        for i in np.arange(len(pole_lons)):
           #Convert central coords to the same as pole coords for consistency
           if verbose: print('Reading pole-list for %s' % (name[i]) )
 
           #If center coords are given, make sure they're in the same system as pole coords
           if (clons[i]!=-999. and clats[i]!=-999.):     
             if  c_coot[i] in pole_coot[i]:
                center=[clons[i],clats[i]]   #If in the same system, just save
                if verbose: print( ' same system, center:',center )
             else:
               if 'gal' in c_coot[i]: clons[i],clats[i]=bovyc.lb_to_radec(clons[i],clats[i],degree=True) 
               else: clons[i],clats[i]=bovyc.radec_to_lb(clons[i],clats[i],degree=True)
               center=[clons[i],clats[i]]    
               if verbose: print(' computed center',center)                     
           else: 
             center=None                          
                                  
           if dlons[i]<-90: dlons[i]=None
           if dlats[i]<0.: dlats[i]=0.                 
                                                       
           #Get realization of great-circle, given stream's orbital pole [and optionally center,length and width]
           gc_lons,gc_lats=gcutils.get_gc_for_pole(pole_lons[i],pole_lats[i],degree=True,center=center,
                                                    dlon=dlons[i],dlat=dlats[i],step=gcstep)
                            
           #Do linear interpolation for the distance, not much better to do
           if ros[i]==rfs[i]: D=ros[i]*np.ones_like(gc_lons)
           elif rfs[i]>ros[i]>0.: D=scipy.interp(gc_lons,[lono[i],lonf[i]],[ros[i],rfs[i]])    
           else: D=-1*np.ones_like(gc_lons)
                                           
           #Create footprint appropriately depending on coord type, this will define the l,b,ra,dec attributes appropriately
           footprint=Footprint(gc_lons,gc_lats,name[i],Dist=D,degree=True,cootype=pole_coot[i])
        
           #Store
           self[name[i]]=footprint     

  #-----------------by lonlat range--------------------------------------------------------------
  def init_by_lonlat_range(self,libname=None,Nran=1000,verbose=False,Rstat='mean'):
                
    if libname is None:
     lib_llrange_filen=os.path.join(os.path.dirname(os.path.realpath(__file__)),'lib','lib_by_lonlat_range.dat')
    else: lib_llrange_filen=libname

    #Validate options for Rstat
    stat_types = ['min', 'max', 'mean', 'median']
    if Rstat not in stat_types:
        raise ValueError("Invalid stat type. Expected one of: %s" % stat_types)
    Rstat_func=getattr(np,Rstat)

    
    name,cootype=scipy.genfromtxt(lib_llrange_filen,usecols=(0,7),unpack=True,dtype=str)
    azo_l,azf_l,lato_l,latf_l,ro_l,rf_l,stype=scipy.genfromtxt(lib_llrange_filen,usecols=(1,2,3,4,5,6,8),unpack=True,
                                                               filling_values=-999.)

    ind=np.arange(azo_l.size)

    for i in ind:
     #Parse corners
     if verbose: print('Reading lon-lat range for %s' % (name[i]) )
     azo,azf,lato,latf,ro,rf=azo_l[i],azf_l[i],lato_l[i],latf_l[i],ro_l[i],rf_l[i]

     #Random realization of coords    
     Rs,azs,lats=get_random_spherical_coords(Nran,rad=[ro,rf],az=[azo,azf],lat=[lato,latf],degree=True)

     #Create footprint appropriately depending on coord type, this will define the l,b,ra,dec attributes appropriately
     footprint=Footprint(azs,lats,name[i],Dist=Rstat_func(Rs)*np.ones_like(Rs),degree=True,cootype=cootype[i])
        
     #Store
     self[name[i]]=footprint  
        

  #-----------------------by list of stars-------------------------------------------------------
  def init_by_star_list(self,verbose=False,Nint=300):
      
     lib_path=os.path.join(os.path.dirname(os.path.realpath(__file__)),'lib')      
     lib_by_stars_list_filen=os.path.join(lib_path,'lib_by_star.log')
            
     list_dat=scipy.genfromtxt(lib_by_stars_list_filen,usecols=(1-1,2-1,3-1,4-1,5-1,6-1))
     list_sdat=scipy.genfromtxt(lib_by_stars_list_filen,usecols=(7-1,8-1,9-1),dtype=str)       
        
     #Deal with one-liners   
     if np.ndim(list_dat)==1: 
       list_dat,list_sdat=np.reshape(list_dat,(1,list_dat.size)),np.reshape(list_sdat,(1,list_sdat.size))
            
     #Parse       
     lon_col,lat_col,Rhel_col,interp_flg,progtype_col,cntflg_col=list_dat.astype(int).T      
     #Convert col numbers to python-style
     lon_col,lat_col,Rhel_col,progtype_col,cntflg_col=lon_col-1,lat_col-1,Rhel_col-1,progtype_col-1,cntflg_col-1       
     cootype,name,fname=list_sdat.T    

     #Loop over library of star list files
     for i in range(list_dat[:,0].size):
       
       if verbose: print('Reading star list for %s' % (name[i]))
       azs,lats,Rhels,cntflg=scipy.genfromtxt(os.path.join(lib_path,fname[i]),usecols=(lon_col[i],lat_col[i],
                                              Rhel_col[i],cntflg_col[i]),unpack=True) 
       
       if interp_flg[i]==1 and azs.size>1:
         #connect these with linear interpolation  
         az_int=np.linspace(azs.min(),azs.max(),Nint)
         lat_spline=scipy.interpolate.CubicSpline(azs,lats,bc_type='natural')
         Rhel_spline=scipy.interpolate.CubicSpline(azs,Rhels,bc_type='natural')
         lat_int,Rhel_int=lat_spline(az_int),Rhel_spline(az_int)

         #If any points have Rhel=-1, make sure not to use them for interpolation
         if (Rhels<0).any():
          #range of azs with valid Rhel:
          azmin,azmax=np.min(azs[Rhels>=0]),np.max(azs[Rhels>=0])
          Rhel_int[(az_int<azmin) | (az_int>azmax)]=-1
         #If there are negative latitudes given, make sure to translate to [0,360] 
         if (az_int<0).any():
          az_int[az_int<0] = 360.+ az_int[az_int<0]    
       else:
        az_int,lat_int,Rhel_int=azs,lats,Rhels          
                          
       #Create footprint appropriately depending on coord type, this will define the l,b,ra,dec attributes appropriately
       footprint=Footprint(az_int,lat_int,name[i],Dist=Rhel_int,degree=True,cootype=cootype[i])
                        
       #If cntflg_col==1, force center to this coords (instead of the default vector mean)
       if cntflg_col[i]>=0:
        if (cntflg==1).sum()>1: print('Warning, more than one point flaged as center in %s, choosing first one.' % (fname[i]))
        caz,clat=azs[cntflg==1][0],lats[cntflg==1][0]    
        auxfoot=Footprint(caz,clat,'dummy',degree=True,cootype=cootype[i]) 
        footprint.cra,footprint.cdec=auxfoot.cra,auxfoot.cdec
        footprint.cl,footprint.cb=auxfoot.cl,auxfoot.cb
        
       #Store
       self[name[i]]=footprint       

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

  def __init__(self,verbose=True,gcstep=0.1,N=1000,Rstat='mean'):
        
    #A MWStreams object is a dictionary in which each entry is a Footprint object, indexed by each stream's name.
    #Since everyone likes to define a stream in a different way, there are several contructor methods to create the stream's
    #footprints, depending on what data is available in the literature
      
    #Initialize empty dictionary
    dict.__init__(self) 
    
    #---For streams defined by RA-DEC or l-b range, create and append footprints
    self.init_by_lonlat_range(verbose=verbose,Nran=N,Rstat=Rstat) 
    
    #---For streams defined by their end-points, create and append footprints
    self.init_by_end_points(verbose=verbose,gcstep=gcstep)
    
    #---For streams defined by their poles, create and append footprints
    self.init_by_poles(verbose=verbose,gcstep=gcstep)   
         
    #---For streams defined by an explicit list of stars, create and append footprints    
    self.init_by_star_list(verbose=verbose)
    
    #---Make sure galactocentric attributes are set to None if Rhel<0
    for i in self.keys():
      #If there's no distance info at all, set to none
      if (self[i].Rhel<0.).all(): 
        self[i].Rgal,self[i].phi,self[i].theta,self[i].cphi,self[i].ctheta=None,None,None,None,None 
      #If some part has no distance info (e.g. the southern extension of the Orphan stream), set that part to null value
      elif (self[i].Rhel<0.).any(): 
        self[i].Rgal[self[i].Rhel<0.]=np.nan
        self[i].phi[self[i].Rhel<0.]=np.nan
        self[i].theta[self[i].Rhel<0.]=np.nan


    #---Override default labelling coords when pre-defined values available.
    self.load_user_defined_centers_and_shortnames()

  #-------------method to plot whole MW streams compilation object at once------------------------------------
  def plot_stream_compilation(self,ax,Rstat='mean',Rrange=[0.,9e9],plot_stream_type='all',plot_names=True,
                              use_shortnames=False,plot_colorbar=False,
                              scat_kwargs=None,text_kwargs=None,sym_kwargs=None,cb_kwargs=None,cootype='gal',verbose=False,
                              exclude_streams=[],include_only=[]): 

   '''plot_stream_compilation(self,ax,Rstat='mean',Rrange=[0.,9e9],plot_stream_type='all',plot_names=True,
                              use_shortnames=False,plot_colorbar=False,
                              scat_kwargs=None,text_kwargs=None,sym_kwargs=None,cb_kwargs=None,cootype='gal',verbose=False,
                              exclude_streams=[],include_only=[])  '''


   #Validate options for cootype
   coo_types = ['gal', 'equ', 'GC']
   if cootype not in coo_types:
        raise ValueError("Invalid coo type. Expected one of: %s" % coo_types)
        
   #Validate options for Rstat
   stat_types = ['min', 'max', 'mean', 'median']
   if Rstat not in stat_types:
        raise ValueError("Invalid stat type. Expected one of: %s" % stat_types)
   Rstat_func=getattr(np,Rstat)

   #Validate options for stream types to be plotted
   stream_types = ['all', 'gc', 'dwarf']
   if plot_stream_type not in stream_types:
        raise ValueError("Invalid stream type. Expected one of: %s" % stream_types)

   #Set a few plotting and labelling defaults  
   Rmax=0.
   if 'GC' in cootype: 
      for i in self.keys(): Rmax=np.max([Rmax,np.max(self[i].Rgal)]) 
   else: 
      for i in self.keys(): Rmax=np.max([Rmax,np.max(self[i].Rhel)]) 
   scatter_kwargs=dict(marker='o',s=8.,edgecolor='none',vmin=0.,vmax=Rmax,alpha=0.5)
   textlab_kwargs=dict(horizontalalignment='center',verticalalignment='bottom',zorder=99)
   textsym_kwargs=dict(marker='+',color='k',ms=5,zorder=textlab_kwargs['zorder'])
   colorbar_kwargs=dict(extend='max')
   if Rrange[1]<9e9: 
     colorbar_kwargs=dict(extend='max')
     if Rrange[0]>0: colorbar_kwargs=dict(extend='both')
   elif Rrange[0]>0: colorbar_kwargs=dict(extend='min')

   #but override whichever are user-supplied (doing it this way I ensure using my own defaults and not matplotlib's
   #if user supplies values for some (but not all) keywords
   if scat_kwargs is not None:
        for key in scat_kwargs.keys(): scatter_kwargs[key]=scat_kwargs[key]
   if text_kwargs is not None:
        for key in text_kwargs.keys(): textlab_kwargs[key]=text_kwargs[key]
   if sym_kwargs is not None:
        for key in sym_kwargs.keys(): textsym_kwargs[key]=sym_kwargs[key]
   if cb_kwargs is not None:
        for key in cb_kwargs.keys(): colorbar_kwargs[key]=cb_kwargs[key]       

   #If Rrange is supplied, use it to set vmax 
   if Rrange[1]<9e9: scatter_kwargs['vmax']=Rrange[1]

   #------------------------------PLOT---------------------------------------------------------------------- 
   keys_to_plot=self.keys()
   if include_only: keys_to_plot=include_only

   for i in keys_to_plot:

    #Skip it stream name in list of excluded streams    
    if i in exclude_streams:
        print('Skipping excluded stream: %s' % (self[i].name))
        continue

    #Skip it if coo-mode selected is galactocentric and stream has no valid galactocentric attributes  
    if 'GC' in cootype and self[i].phi is None:
        print('Skipping stream %s, no valid Rhel => no valid galactocentric attributes' % (self[i].name))
        continue
       
    if 'GC' in cootype:
      notnan=self[i].Rgal>=0. #this works to avoid nans as well
      ro,rf,Rs=np.min(self[i].Rgal[notnan]),np.max(self[i].Rgal[notnan]),self[i].Rgal    
    else:
      notneg=self[i].Rhel>=0.
      if notneg.any():
        ro,rf,Rs=np.min(self[i].Rhel[notneg]),np.max(self[i].Rhel[notneg]),self[i].Rhel    
      else: ro,rf,Rs=-1,-1,self[i].Rhel

    if 'gal' in cootype: lons,latts,loncenter,latcenter=self[i].l,self[i].b,self[i].cl,self[i].cb
    elif 'equ' in cootype: lons,latts,loncenter,latcenter=self[i].ra,self[i].dec,self[i].cra,self[i].cdec      
    else: lons,latts,loncenter,latcenter=self[i].phi,self[i].theta,self[i].cphi,self[i].ctheta 
    

    #Skip if stream does not have any part overlapping Rrange
    if (ro>0) and not (ro<=Rrange[1] and rf>=Rrange[0]):
      if verbose: 
        print('Skipping %s [%.1f,%.1f], outside selected Rrange [%.1f,%.1f]' % (self[i].name,ro,rf,Rrange[0],Rrange[1]))
      continue
    else:
      if verbose and ro is not None:
        print('Plotting %s [%.1f,%.1f], INside selected Rrange [%.1f,%.1f]' % (self[i].name,ro,rf,Rrange[0],Rrange[1]))
   
    cc=ax.scatter(lons,latts,c=Rs,**scatter_kwargs)
    if plot_names:
     if use_shortnames: plotname=self[i].sname
     else:              plotname=self[i].name
     ax.plot(loncenter,latcenter,**textsym_kwargs)
     ax.text(loncenter,latcenter,plotname,**textlab_kwargs)
     
    if cc and plot_colorbar: 
        cbar=plt.colorbar(cc,ax=ax,**colorbar_kwargs)
        cbar.solids.set_rasterized(True)   #To avoid white lines in colorbar that appear some times
        cbar.solids.set_edgecolor("face")
        plot_colorbar=False
  

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
