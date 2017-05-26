import numpy as np
import scipy
import pylab as plt
import mypyutils as mylib
import myutils
import bovy_coords as bovyc
import gcutils 

def plot_PS1_streams(ax,plot_colorbar=False,scat_kwargs=None,text_kwargs=None,sym_kwargs=None,galactic=True,gcstep=0.1):
    
 #=======PS1 stream data from Bernard+ 2016(arXiv:1607.06088)=====================================================
 ps_labels=['PS1-A','PS1-B','PS1-C','PS1-D','PS1-E']
 ps_pole_ras=np.array([300.856,65.603,232.227,49.640,42.526])
 ps_pole_decs=np.array([20.732,32.567,33.838,2.467,23.987])
 ps_clons=[160.17,248.41,75.12,231.06,144.17] #stream's centre coords given in galactic
 ps_clats=[-62.27,32.30,-32.6,32.83,58.40]
 ps_lens=[5.,10.,8.,45.,25.]
 ps_widths=np.array([27,27,20,52,37])/60.
 ps_rhels=[7.9,14.5,17.4,22.9,12.6]

 #set a few plotting and labelling defaults  
 scatter_kwargs=dict(s=10,vmin=0.,edgecolor='none')
 textlab_kwargs=dict(horizontalalignment='center',verticalalignment='bottom',zorder=99)
 textsym_kwargs=dict(marker='+',color='k',ms=5,zorder=textlab_kwargs['zorder'])      

 #but override whichever are user-supplied (doing it this way I ensure using my own defaults and not matplotlib's
 #if user supplies values for some (but not all) keywords
 if scat_kwargs is not None: 
        for key in scat_kwargs.keys(): scatter_kwargs[key]=scat_kwargs[key]  
 if text_kwargs is not None: 
        for key in text_kwargs.keys(): textlab_kwargs[key]=text_kwargs[key]           
 if sym_kwargs is not None:
        for key in sym_kwargs.keys(): textsym_kwargs[key]=sym_kwargs[key]           
    
 #Loop over PS-1 streams    
 for ii in np.arange(len(ps_pole_ras)):
  ps_pole_ra,ps_pole_dec=ps_pole_ras[ii],ps_pole_decs[ii]
  #Convert central coords to equatorial
  cra,cdec=bovyc.lb_to_radec(ps_clons[ii],ps_clats[ii],degree=True,epoch=2000.0)   
  #Get random realization of great-circle, given stream's orbital pole
  ps_ra,ps_dec=gcutils.get_gc_for_pole(ps_pole_ra,ps_pole_dec,degree=True,center=[cra,cdec],dlon=ps_lens[ii],
                               dlat=ps_widths[ii],step=gcstep)   

  #convert to galactic    
  if galactic: ps_l,ps_b=bovyc.radec_to_lb(ps_ra,ps_dec,degree=True).T
  else: 
    ps_l,ps_b=ps_ra,ps_dec    
    ps_clons[ii],ps_clats[ii]=cra,cdec
    
  #Add layer to plot     
  cc=ax.scatter(ps_l,ps_b,marker='o',c=ps_rhels[ii]*np.ones_like(ps_l),**scatter_kwargs)
  ax.text(ps_clons[ii],ps_clats[ii],ps_labels[ii],**textlab_kwargs)
  ax.plot([ps_clons[ii],],[ps_clats[ii],],**textsym_kwargs)      

  if plot_colorbar: plt.colorbar(cc,ax=ax)          
            
 return cc      

#--------------------------------

def plot_orphan_stream(ax,plot_colorbar=False,scat_kwargs=None,text_kwargs=None,galactic=True):
  
  #set a few plotting and labelling defaults  
  scatter_kwargs=dict(marker='H',s=80,vmin=0.,edgecolor='none',zorder=15)
  textlab_kwargs=dict(horizontalalignment='center',fontsize=13,verticalalignment='bottom',zorder=99)

  #but override whichever are user-supplied (doing it this way I ensure using my own defaults and not matplotlib's
  #if user supplies values for some (but not all) keywords
  if scat_kwargs is not None: 
        for key in scat_kwargs.keys(): scatter_kwargs[key]=scat_kwargs[key]  
  if text_kwargs is not None: 
        for key in text_kwargs.keys(): textlab_kwargs[key]=text_kwargs[key]           
    
  scatter_kwargs2=scatter_kwargs.copy()  
  scatter_kwargs2['s']=scatter_kwargs['s']/8.

  #Newberg 2010 detections for Orphan stream
  orphan_l,orphan_b,orphan_R=scipy.genfromtxt('newberg2010.orphan_orbit.dat',unpack=True)    
  if not galactic: orphan_l,orphan_b=bovyc.lb_to_radec(orphan_l,orphan_b,degree=True).T  
    
  #Plot Orphan stream layer
  #plot the actual measurements  
  ax.scatter(orphan_l,orphan_b,c=orphan_R,**scatter_kwargs)    
  #connect these with linear interpolation  
  orphan_l_int=np.linspace(orphan_l.min(),orphan_l.max(),100)
  orphan_b_int=scipy.interp(orphan_l_int,orphan_l,orphan_b)
  orphan_R_int=scipy.interp(orphan_l_int,orphan_l,orphan_R)
  cc=ax.scatter(orphan_l_int,orphan_b_int,c=orphan_R_int,**scatter_kwargs2)
                #marker='o',s=10,vmin=0.,vmax=Rmax,edgecolor='none',cmap=my_cmap,zorder=15)    

  if plot_colorbar: plt.colorbar(cc,ax=ax)                
        
  #plot Orphan label  
  ax.text(orphan_l[-2],orphan_b[-2],'Orphan',**textlab_kwargs) 
    
  #Orphan stream extension reported by GrillmairXXX  
  orphan_ext_delta=np.linspace(-18.,-38.,100)
  orphan_ext_alpha=163.147-0.0896*orphan_ext_delta+0.00804*(orphan_ext_delta**2)
  orphan_ext_l,orphan_ext_b=bovyc.radec_to_lb(orphan_ext_alpha,orphan_ext_delta,degree=True).T
  if galactic: ax.plot(orphan_ext_l,orphan_ext_b,'-',color='gray')
  else: ax.plot(orphan_ext_alpha,orphan_ext_delta,'-',color='gray')
    

#-----------------------------------------

def plot_globular_clusters(ax,plot_colorbar=False,scat_kwargs=None,galactic=True):
    
 #Harris's Globular cluster compilation
 gc_l,gc_b,gc_Rhel=scipy.genfromtxt('/Users/cmateu/trabajo/local_group/globular_clusters_compilations/globular_cluster_params.harris.tableI.csv',
                              missing_values="",usecols=(5-1,6-1,7-1),unpack=True,delimiter=',',dtype=np.float)
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
       
#------------------------------------------

def plot_grillmair_carlin_streams(ax,Rstat='median',scat_kwargs=None,text_kwargs=None,sym_kwargs=None,
                                  galactic=True,plot_stream_type='all',plot_names=True,Nran=2000):

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
 scatter_kwargs=dict(marker='o',s=8.,edgecolor='none',vmin=0.,alpha=0.5)
 textlab_kwargs=dict(horizontalalignment='center',verticalalignment='bottom',zorder=99)
 textsym_kwargs=dict(marker='+',color='k',ms=5,zorder=textlab_kwargs['zorder'])      

 #but override whichever are user-supplied (doing it this way I ensure using my own defaults and not matplotlib's
 #if user supplies values for some (but not all) keywords
 if scat_kwargs is not None: 
        for key in scat_kwargs.keys(): scatter_kwargs[key]=scat_kwargs[key]  
 if text_kwargs is not None: 
        for key in text_kwargs.keys(): textlab_kwargs[key]=text_kwargs[key]           
 if sym_kwargs is not None:
        for key in sym_kwargs.keys(): textsym_kwargs[key]=sym_kwargs[key]           
   
 #==========================Known Streams Grillmair & Carlin===========================================================
 azo_l,azf_l,lato_l,latf_l,ro_l,rf_l,stype=scipy.genfromtxt('grillmair_carlin2016_known_streams_full.dat',usecols=(1,2,3,4,5,6,8),unpack=True)
 name,cootype=scipy.genfromtxt('grillmair_carlin2016_known_streams_full.dat',dtype='S',usecols=(0,8-1),unpack=True)

 ind=np.arange(azo_l.size)
   
 for i in ind:
              
   #Skip type depending on selected options
   if 'gc' in plot_stream_type and stype[i]==2.: continue     #If on 'gc' mode, skip dwarf galaxy streams
   if 'dwarf' in plot_stream_type and stype[i]==1.: continue  #If on 'dwarf' mode, skip globular cluster streams               
                
   #Parse corners   
   azo,azf,lato,latf,ro,rf=azo_l[i],azf_l[i],lato_l[i],latf_l[i],ro_l[i],rf_l[i]    
    
   #Random realization of coords-polygon    
   Rs,azs,lats=mylib.get_random_spherical_coords(Nran,rad=[ro,rf],az=[azo,azf],lat=[lato,latf],degree=True)  
   #Center coords   
   az_center,lat_center,R_selected=np.mean(azs),np.mean(lats),Rstat_func(Rs)   
        
   
   if 'gal' in cootype[i] and galactic:       #If coords given are galactic and plotting mode is galactic
     lons,lats=azs,lats
     loncenter,latcenter=az_center,lat_center
   elif 'gal' in cootype[i] and not galactic: #If coords given are galactic and plotting mode is equatorial
     lons,lats=bovyc.lb_to_radec(azs,lats,degree=True).T 
     loncenter,latcenter=bovyc.lb_to_radec(az_center,lat_center,degree=True)           
   elif 'eq' in cootype[i] and galactic:     #If coords given are equatorial and plotting mode is galactic
     lons,lats=bovyc.radec_to_lb(azs,lats,degree=True).T 
     loncenter,latcenter=bovyc.radec_to_lb(az_center,lat_center,degree=True)           
   #else:                                  #If coords given are equatorial and plotting mode is equatorial
   elif 'eq' in cootype[i] and not galactic:
     lons,lats=azs,lats
     loncenter,latcenter=az_center,lat_center           
 
   #------------------------------PLOT----------------------------------------------------------------------  
   ax.scatter(lons,lats,c=R_selected*np.ones_like(lons),**scatter_kwargs)  
   if plot_names:
    ax.plot(loncenter,latcenter,**textsym_kwargs)
    ax.text(loncenter,latcenter,name[i],**textlab_kwargs)   

#----------------------------------------------



