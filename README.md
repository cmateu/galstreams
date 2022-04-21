# **galstreams**

![see plot here](examples/fig_all_streams_lib.png?raw=true "galstreams 04-2022")

### DESCRIPTION:

The new and improved *galstreams* Library of Stellar Streams in the Milky Way (v1.0) is introduced. Stellar streams are now supported as track SkyCoord objects (Track6D), rather than the Footprint objects provided in the previoues version (v0.1). The main new features of the library are:

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

and source your .cshrc / .bashrc or equivalent file.

If you do not have root access, you can install in the custom directory path_to_dir.
First, add the directory's path path_to_dir and path_to_dir/lib/python??/site-packages/
to the PYTHONPATH variable in your .cshrc (or .bashrc) file and source it. Then install using the --prefix option::

    python setup.py install --prefix=path_to_dir

Add path_to_dir/bin to your PATH in your .csrhc or .bashrc file.

----------
# Quick Guide

A MWstreams object can be easily created as follows:

	mwsts=galstreams.MWStreams()

Running this in verbose mode (verbose=True) will print the each of the library’s stream names as they are initialized (it's getting long, so, perhaps don't). The resulting object, mwsts, contains the footprint information for each of the stream tracks in the library. By default the MWStreams object will only create one default track for each of the distinct stellar streams in the library (see Mateu 2022). To explore all tracks available you can set implement_Off=True.

Each stream track is referenced by it's unique TrackName. Since the MWStreams is just a dictionary, you can list the available TrackNames by doing:

    mwsts.keys()

You can also query the MWStreams object using the get_track_names_for_stream method to find the TrackNames for all the streams that match a given string. For example:

    mwsts.get_track_names_for_stream('Pal')

this will print a list of all the streams that contain the string 'Pal'.

The MWStreams object also has a *summary* attribute, a Pandas DataFrame summarising properties for all the stream tracks in the library:

    mwsts.summary


To make a quick plot of the stream’s library stored in the mwsts object use:

    fig = plt.figure(1,figsize=(16,11))
    ax = fig.add_subplot(111, projection='mollweide')

    for st in mwsts.keys():
      #Plot the tracks  
      ax.scatter(mwsts[st].track.galactic.l.wrap_at(180*u.deg).rad,
                 mwsts[st].track.galactic.b.rad, marker='.', s=30,
                 label="{ID:.0f}={Name}".format(ID=mwsts[st].ID,Name=mwsts.summary.Name[st]))
      #Annotate at one of the end points  
      ax.annotate(mwsts[st].ID, xy=(mwsts[st].end_points.galactic.l.wrap_at(180*u.deg)[0].rad,mwsts[st].end_points.galactic.b[0].rad),
                  xycoords='data',
                  arrowprops=dict(arrowstyle="-",color='k'),
                  horizontalalignment='center', verticalalignment='center',
                  xytext=(-10,15),textcoords='offset points',
                  )

    ax.legend(ncol=8,loc='center', columnspacing=0.5, handletextpad=0.1,
          bbox_to_anchor=(0.5,-0.28), markerscale=3, fontsize='medium')

![see plot here](examples/quickex.png?raw=true "Example plot for galstreams") -->

See the notebooks provided in the examples folder for more detailed examples of the library's functionality.

