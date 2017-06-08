mwstreams
=========

**DESCRIPTION:**

mwstreams is a Python Library and Toolkit for visualizing streams and overdensities known in the Milky Way (MW).

The main contributions of this toolkit are the following:

- The current distribution provides tools to handle and plot each stream's footprint in different coordinate systems. 

- Utility methods are provided to automatically add all MW known globular clusters and dwarf galaxy satellites.

- The MWStreams class provides an object to handle the data for all known streams in the MW

- The gcutils library provides utility methods that will let you create a footprint object for your own stream, defining it in one of the following four ways:

  - by giving the coordinates of the start and end point
  - by giving the orbital pole, and if known its center, length and width 
  - by giving a range of RA/DEC or l/b  (more suitable for clouds rather than for streams)
  - by giving a list of coordinates 

----------

**REQUIREMENTS**

- Python modules required are NUMPY and SCIPY. MATPLOTLIB is needed for plotting utilities.
- This toolkit makes use of the coordinate transformation library bovy_coords.py from the galpy package by Jo Bovy (2015, in prep.). It is supplied with this bundle.

----------

**INSTALLATION**

In a terminal, run the following command::

    sudo python setup.py install 

Source your .cshrc

If you do not have root access, you can install in the custom directory path_to_dir.
First, add the directory's path path_to_dir and path_to_dir/lib/python2.7/site-packages/
to the PYTHONPATH variable in your .cshrc/.bashrc file and source it. Then install using the --prefix option::

    python setup.py install --prefix=path_to_dir

Add path_to_dir/bin to your PATH in your .csrhc or .bashrc file.

*This is not set up yet, but will be soon.*

----------

Quick Guide
===========

A MWstreams object can be easily created as follows::

	mwsts=mwstreams.MWStreams(verbose=True)

Running this in verbose mode will print the each of the library’s stream names as they are initialized.

To quickly plot the the stream’s library stored in it use:: 

	fig=plt.figure(1,figsize=(16,8))
	ax=fig.add_subplot(111)
	mwsts.plot_stream_compilation(ax,plot_colorbar=True)

(plot not shown)

the plot is made in equatorial coordinates by default, but heliocentric galactic and galactocentric spherical coordinates can also be used. For more details on available MWStreams methods see XXX-link-here.

To overplot the positions of MW globular clusters use::

	mwstreams.plot_globular_clusters(ax)

----------

The Built-in Streams Library
============================



----------

Package: mwstreams 
==================

**MWstreams Class**

The MWstreams class returns a dictionary containing a Footprint object (see XXX-link-here) for each of the streams known in the MW, indexed by the stream’s name.

The class can be easily instantiated as follows::

	mwsts=mwstreams.MWStreams()

This will read the stream definitions stored in the lib directory to instantiate a footprint object appropriately for each stream.

**Footprint Class**

This class handles a stream’s footprint as a collection of points. A Footprint object can be instantiated by passing it a name string and a pair of longitude-latitude vectors which will be interpreted as RA/DEC or l/b if the coordinate system is indicated as equatorial or galactic respectively (via de cootype keyword), e.g., as follows::

	footprint= mwstreams.Footprint(ra,dec,’Amethyst’,cootype=‘equ’)

The heliocentric distance, proper motions and radial velocity can be passed at initialization as optional arguments.

The coordinates in other reference systems of interest are computed by default and set as object attributes. 

As an instance of the Footprint class, a footprint object has the following default attributes:

- footprint.ra, footprint.dec    (equatorial coords)
- footprint.l, footprint.b       (galactic coords) 
- footprint.cra, .cdec, .cl, .cb (footprint’s geometric center coordinates)

If Rhel, the heliocentric distance, is given:

- footprint.Rhel
- footprint.phi, .theta       (galactocentric coords)
- footprint.Rgal              (galactocentric distance)
- footprint.xhel,.yhel,.zhel  (cartesian heliocentric coords)
- footprint.x,.y,.z           (cartesian galactocentric coords)

If proper motions are given:

- footprint.pmra, .pmdec, .pmrastar  (pmrastar=pmra*cos(dec))
- footprint.pml, .pmb, .pmlstar       (pmlstar=pml*cos(dec))

If radial velocity is given:

- footprint.vrad

If all above given:

- footprint.vxhel,.vyhel,.vzhel  (cartesian heliocentric vels)
- footprint.vx,.vy,.vz           (cartesian galactocentric vels)

An utility method is provided to apply a mask to all array attributes of a Footprint object::

	Footprint.mask_footprint(mask)

For full details see the doc-string for the Footprint class.

----------

**FILES PROVIDED**

- Executable programs:
	* work in progress - stand-alone code to make a quick plot of the MW library in user-selected coords

- Libraries:
	* mwstreams.py
	* gcutils.py
	* bovy_coords.py
	* pyutils.py
	* lib
- Documentation
   * README.rst








