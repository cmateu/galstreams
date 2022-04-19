# **galstreams**

![see plot here](examples/fig_all_streams_lib.png?raw=true "galstreams 04-2022")

### DESCRIPTION:

A new and improved *galstreams* Library of Stellar Streams in the Milky Way is introduced. Stellar streams are now supported as track SkyCoord objects (Track6D), rather than the Footprint objects provided in v0.1. The main features of the new library are:

-  Celestial, distance, proper motion and radial velocity tracks for each stream (pm/vrad when available) stored as astropy SkyCoord objects
-  Stream's (heliocentric) coordinate frame realised as astropy reference frame
-  Stream's end-points and mid-point
- Polygon Footprints
-  Pole (at mid point) and pole tracks in the heliocentric and Galactocentric (GSR) frames
-  Angular momentum track in a heliocentric reference frame at rest with respect to the Galactic centre
-  Summary object for the full library: Uniformly reported stream length, end points and mid-point, heliocentric and Galactocentric mid-pole, track and discovery references and information flag denoting which of the 6D attributes (sky, distance, proper motions and radial velocity) are available in the track object.

The new library includes 125 stream tracks corresponding to 97 distinct stellar streams (updated as of March 2022). The procedure use for the track  described in detail Mateu 2022.

### REQUIREMENTS

- Python modules required are NUMPY, SCIPY, MATPLOTLIB, ASTROPY and GALA.

----------

### INSTALLATION

In a terminal, run the following command:

    sudo python setup.py install

and source your .cshrc

If you do not have root access, you can install in the custom directory path_to_dir.
First, add the directory's path path_to_dir and path_to_dir/lib/python??/site-packages/
to the PYTHONPATH variable in your .cshrc (or .bashrc) file and source it. Then install using the --prefix option::

    python setup.py install --prefix=path_to_dir

Add path_to_dir/bin to your PATH in your .csrhc or .bashrc file.

----------
# Quick Guide

A MWstreams object can be easily created as follows:

	mwsts=galstreams.MWStreams(verbose=True)

Running this in verbose mode will print the each of the library’s stream names as they are initialized (it's getting long, so, perhaps don't). The resulting object, mwsts, contains the footprint information for each of the stream tracks in the library. By default the MWStreams object will only create the default track for each of the distinct stellar streams in the library (see Mateu 2022). To explore all tracks available you can set implement_Off=True.

Each stream track is referenced by it's unique TrackName. Since the MWStreams is just a dictionary, you can list the available TrackNames by doing:

  mwsts.keys()

The MWStreams object also has a *summary* attribute, a Pandas DataFrame summarising properties for all the stream tracks in the library:

  mwsts.summary.head()



To make a quick plot of the stream’s library stored in the mwsts object use:

	fig=plt.figure(1,figsize=(16,8))
	ax=fig.add_subplot(111)
	cmapp=plt.cm.plasma_r
	cmapp.set_under('grey')   #If distance info is missing (-1), plot in grey
  mwsts.plot_stream_compilation(ax,plot_colorbar=True,scat_kwargs=dict(vmin=0.,vmax=80.,cmap=cmapp, alpha=0.3),
                                      use_shortnames=False, cb_kwargs=dict(label='Heliocentric Distance (kpc)'),
                                      verbose=False)


	ax.set_xlim(0.,360.)
	ax.set_ylim(-90.,90.)
	ax.set_xlabel('$l$ (deg)')
	ax.set_ylabel('$b$ (deg)')

![see plot here](examples/quickex.png?raw=true "Example plot for galstreams")

The plot is made in galactic coordinates by default, but equatorial and galactocentric spherical coordinates can also be used (cootype="equ" or cootype="gc") ![see plots here](examples/quickex_ra_phitheta.png?raw=true). This example shows how you can make changes to the plot by passing arguments to the 'scatter' and 'colorbar' commands through the scat_kwargs and cb_kwargs keyword arguments. Also, axis text and symbol-text (stream names) properties can be modified using the text_kwargs and sym_kwargs. For more details on available MWStreams methods see [here](#mwstreams-class).

MW globular clusters data from the Harris (1996, 2010 edition) compilation is also included in the library. To quickly overplot globular clusters as an extra layer in the previous add:

	galstreams.plot_globular_clusters(ax)

There are several options available to customize these plots, to check them out have a look at the MWStreams doc-string.

Here's also an example on how to use the stream's reference frame attribute. The plot below shows several streams that intersect the GD-1 track, all plotted in GD-1's reference frame. This can be produced with the following code:

        #import seaborn as sns    #---uncomment these two lines to match the style of the plot shown below---
        #sns.set_context('talk')
	plt.figure(1,figsize=(12,5))
	plt.subplot(111)
	for ss in ['GD-1','Gaia-5','PS1-D','PS1-E','Orphan']:
         #The reference frame attribute (gcfr) for GD-1 is used to convert each stream's astropy.SkyCoord object
         #to GD-1's reference frame
	 scoo = mwsts[ss].sc.transform_to(mwsts['GD-1'].gcfr)   
	 plt.plot(scoo.phi1,scoo.phi2,'.',label=ss)

	plt.legend(ncol=2,markerscale=2)
	plt.xlim(-60,60)
	plt.ylim(-20,20)
	plt.xlabel('$\phi_1$ ($\degree$)')
	plt.ylabel('$\phi_2$ ($\degree$)')
	plt.title("Streams in GD-1's reference system" )
	plt.tight_layout()

![see plot here](examples/quickex_gd1_ref_system.png?raw=true "Example plot for galstreams")
