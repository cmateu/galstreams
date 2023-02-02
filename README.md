# **galstreams**

![see plot here](notebooks/fig_all_streams_lib.png?raw=true "galstreams 03-2022")

### DESCRIPTION:

The new and improved *galstreams* Library of Stellar Streams in the Milky Way (v1.0) is introduced. Stellar streams are now supported as track SkyCoord objects (Track6D), rather than the Footprint objects provided in the previoues version (v0.1). The main new features of the library are:

-  Celestial, distance, proper motion and radial velocity tracks for each stream (pm/vrad when available) stored as astropy SkyCoord objects
-  Stream's (heliocentric) coordinate frame realised as astropy reference frame
-  Stream's end-points and mid-point
- Polygon Footprints
-  Pole (at mid point) and pole tracks in the heliocentric and Galactocentric (GSR) frames
-  Angular momentum track in a heliocentric reference frame at rest with respect to the Galactic centre
-  Summary object for the full library: Uniformly reported stream length, end points and mid-point, heliocentric and Galactocentric mid-pole, track and discovery references and information flag denoting which of the 6D attributes (sky, distance, proper motions and radial velocity) are available in the track object.

The new library includes 126 stream tracks corresponding to 95 distinct stellar streams (updated as of March 2022). The library is described in detail in [Mateu (2023)](https://arxiv.org/abs/2204.10326).

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

For a detailed walk through the library please see the example Python notebooks provided [here](notebooks/).
