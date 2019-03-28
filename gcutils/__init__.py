import numpy as np
import bovy_coords as bovyc
import galstreams

def great_circle_distance(lono,lato,lonf,latf,degree=False):

    if degree: f=np.pi/180.
    else: f=1.

    return (1./f)*np.arccos(np.sin(latf*f)*np.sin(lato*f) + np.cos(latf*f)*np.cos(lato*f)*np.cos(f*(lonf-lono))) 

def get_gc_for_pole(_lon,_lat,degree=True,step=0.01,dlat=0.,center=None,dlon=None,Nran=10000):
    if degree: f=np.pi/180.
    else: f=1.
    lon,lat=f*_lon,f*_lat

    #Generate great circle with pole = u_z
    azs,thetas=np.array([]),np.array([])
#    n = np.int(Nran/step)
#    azs,thetas=galstreams.get_random_spherical_angles(n,az=[0.,2*np.pi],lat=[-f*dlat/2.,+f*dlat/2.],degree=False) 
    if dlat>0:
     for lato in np.radians(np.arange(-dlat/2.,dlat/2.,step)):
        aux=np.radians(np.arange(0.,360.,step))
        azs=np.append(azs,aux)
        thetas=np.append(thetas,lato*np.ones_like(aux))
    else:    
     azs=np.radians(np.arange(0.,360.,step))
     thetas=np.zeros_like(azs) 
        
    _x=np.cos(azs)*np.cos(thetas)
    _y=np.sin(azs)*np.cos(thetas)
    _z=np.sin(thetas)     
    
    #Rotate so that pole has now the given lon,lat
    coslon,sinlon,coscolat,sincolat=np.cos(lon),np.sin(lon),np.cos(np.pi/2.-lat),np.sin(np.pi/2.-lat)
    _x1=  coscolat*_x + sincolat*_z
    _y1= _y
    _z1= -sincolat*_x + coscolat*_z
    
    _x2= coslon*_x1 - sinlon*_y1
    _y2= sinlon*_x1 + coslon*_y1
    _z2=_z1
    
    m=bovyc.XYZ_to_lbd(_x2,_y2,_z2,degree=degree)
    phi2,theta2,Rgal2=m.T

    if center and dlon:
        cphi,ctheta=center
        dist=(1./f)*np.arccos(np.sin(theta2*f)*np.sin(ctheta*f) + np.cos(theta2*f)*np.cos(ctheta*f)*np.cos(f*(phi2-cphi)))     
        mask=dist<=dlon/2.   
        phi2,theta2=phi2[mask],theta2[mask]
    
    return (phi2,theta2)


def get_gc_for_pair(lono,lato,lonf,latf,step=0.01,degree=True,dlat=0.):

   #Find pole cartesian coords
   u_ini=bovyc.lbd_to_XYZ(lono,lato,1.,degree=degree)
   u_end=bovyc.lbd_to_XYZ(lonf,latf,1.,degree=degree)
   u_ini,u_end=np.array(u_ini),np.array(u_end)

   #Do cross product to find great-circle pole
   pX,pY,pZ=np.cross(u_ini,u_end)

   #Convert pole coords back to spherical
   pole_lon,pole_lat=bovyc.XYZ_to_lbd(pX,pY,pZ,degree=degree)[:2]

   #Find center as average of the two vectors
   cX,cY,cZ=u_ini+u_end   #these are cartesian, so we're good
   clon,clat=bovyc.XYZ_to_lbd(cX,cY,cZ,degree=degree)[:2]
  
   #Length equals great circle distance between end-points
   dlon=great_circle_distance(lono,lato,lonf,latf,degree=degree)
 
   #Get great circle arc connecting the given end-points
   gc_lons,gc_lats=get_gc_for_pole(pole_lon,pole_lat,degree=degree,step=step,dlat=dlat,dlon=dlon,center=[clon,clat])

   return gc_lons,gc_lats,clon,clat

