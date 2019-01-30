# galstreams


### DESCRIPTION:

galstreams is a Milky Way Streams Footprint Library and Toolkit for Python. It includes a compilation of spatial information for known stellar streams and overdensities in the Milky Way (MW) and Python tools for visualizing them. It is introduced in [Mateu, Read & Kawata (2018)](http://adsabs.harvard.edu/abs/2018MNRAS.474.4112M).

The main contributions of this package are the following:

- It provides tools to easily handle and plot a stream's footprint in different celestial coordinate systems.

- The MWStreams class provides an object to handle the data for all known streams in the MW at once

- It is easily extendable to include new streams as they are published

- The gcutils library provides utility methods that will let you create a footprint object for your own stream, using any of the four available constructors:

  - by giving the coordinates of the start and end point
  - by giving the orbital pole and, if known, its center, length and width
  - by giving a range of RA,DEC or l,b  (more suitable for clouds rather than for streams)
  - by giving a list of coordinates (RA,DEC or l,b)

- Utility methods are also provided to automatically add all MW known globular clusters and dwarf galaxy satellites.

- The streams_lib_notes.ipynb iPython Notebook keeps record of intermediate computations made for some of the stream's data to comply with one the four available constructors.

The galstreams package can be used to replicate Figures 2, 3 and 4 in [Mateu, Read & Kawata (2018).](http://adsabs.harvard.edu/abs/2018MNRAS.474.4112M)

----------

### VERSION HISTORY:

- 2019/01/30: Phlegethon stream added (Ibata et al. 2019)
- 2019/01/30: Orphan stream track updated (Koposov et al. 2019)
- 2019/01/30: Gaia1-5 streams (Malhan et al. 2018) added
- 2018/04/24: Tucana III stream (Drlica-Wagner et al. 2015, Shipp et al. 2018) added
- 2018/04/23: Compilation of MW satellites added
- 2018/01/24: Shortnames changed for Turbio (Trb) and Turramburra (Trn)
- 2018/01/19: Corvus and 20.0-1 stream's (Mateu et al. 2018 high-confidence candidates) footprints added 
- 2018/01/19: New functionality added: Footprints can be initialized with cootype=GC 
- 2018/01/19: Orinoco and Murrumbidgee footprints corrected. Shortname attribute added. Use shortname option in plots added.
- 2018/01/17: Label centering changed. Now can be set arbitrarily at log file. Labels checked in gal/equ/GC
- 2018/01/17: *New DES streams* (Shipp et al. 2018) *and Jet stream* (Jethwa et al. 2017) *added*.
	      Width option added to end-points constructor and now included for all streams defined this way. 
- 2017/11/21: Eri/Phe coordinates fixed and reference added to Table.

### REQUIREMENTS

- Python modules required are NUMPY, SCIPY and MATPLOTLIB.
- This toolkit makes use of the coordinate transformation library bovy_coords.py from the [galpy](http://github.com/jobovy/galpy) package by [Jo Bovy (2015)](http://adsabs.harvard.edu/abs/2015ApJS..216...29B). It is supplied with this bundle.

----------

### INSTALLATION

In a terminal, run the following command:

    sudo python setup.py install

and source your .cshrc

If you do not have root access, you can install in the custom directory path_to_dir.
First, add the directory's path path_to_dir and path_to_dir/lib/python2.7/site-packages/
to the PYTHONPATH variable in your .cshrc (or .bashrc) file and source it. Then install using the --prefix option::

    python setup.py install --prefix=path_to_dir

Add path_to_dir/bin to your PATH in your .csrhc or .bashrc file.

----------
# Quick Guide

A MWstreams object can be easily created as follows:

	mwsts=galstreams.MWStreams(verbose=True)

Running this in verbose mode will print the each of the library’s stream names as they are initialized. The resulting object, mwsts, contains the footprint information for each of the streams and clouds registered in the library.

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

----------

# The Milky Way Streams Library

The data for the built-in streams library is stored in the *lib* directory, where four main files are found, depending on which constructor is used in galstreams.MWStreams to
define a stream’s footprint:

- lib_by_pair.dat (input: coordinates of the start and end point)
- lib_by_pole.dat (input: orbital pole, and optionally center, length and width)
- lib_by_lonlat_range.dat (input: range of RA/DEC or l/b)
- lib_by_star.log (input: list of star coordinates. this is a log file where the format and location of the star list files are set for each defined stream)  


The following table summarizes the streams included in the library. The list of streams is based in the Grillmair & Carlin (2016) review (their Table 4.1) and updated as of 17/Jan/2018.


| Name         | Reference                  |   Name      |  Reference                       |
|--------------|----------------------------|-------------|----------------------------------|
| Alpheus      |  Grillmair 2013            |   PAndAS    |  Grillmair & Carlin 2016         |
| Acheron      |  Grillmair 2009            |   Phoenix   |  Balbinot 2016                   |
| ACS          |  Grillmair 2006            |   PiscesOv  |  Grillmair & Carlin 2016         |
| ATLAS        |  Koposov 2014              |   PS1-A     |  Bernard 2016                    |
| Cetus        |  Newberg 2009              |   PS1-B     |  Bernard 2016                    |
| Cocytos      |  Grillmair 2009            |   PS1-C     |  Bernard 2016                    |
| GD-1         |  Grillmair 2006            |   PS1-D     |  Bernard 2016                    |
| EBS          |  Grillmair & Carlin 2016   |   PS1-E     |  Bernard 2016                    |
| Eridanus     |  Myeong2017                |   Sagitarius|  Law & Majewski 2010 (model)     |
| Eri/Phe      |  Li2016                    |   Sangarius |  Grillmair 2017                  |
| Hermus       |  Grillmair 2014            |   Scamander |  Grillmair 2017                  |
| Her-Aq       |  Grillmair & Carlin 2016   |   Styx      |  Grillmair 2009                  |
| Hyllus       |  Grillmair 2014            |   Tri-And   |  Grillmair & Carlin 2016         |
| Kwando       |  Grillmair 2017b           |   Tri-And2  |  Grillmair & Carlin 2016         |
| Lethe        |  Grillmair 2009            |   Tri/Pis   |  Bonaca 2012                     |
| Molonglo     |  Grillmair 2017b           |   VOD/VSS   |  Grillmair & Carlin 2016         |
| Monoceros    |  Grillmair & Carlin 2016   |   WG1       |  Agnello 2017                    |
| Murrumbidgee |  Grillmair 2017b           |   WG2       |  Agnello 2017                    |
| NGC5466      |  Grillmair & Johnson2006   |   WG3       |  Agnello 2017                    |
| Ophiucus     |  Bernard 2014              |   WG4       |  Agnello 2017                    |
| Orphan       |  Newberg 2010              |   Jet       |  Jethwa et al. 2017              |
| Orinoco      |  Grillmair 2017b           |   Indus     |  Shipp et al. 2018               |
| Pal5         |  Grillmair 2006            |   Jhelum    |  Shipp et al. 2018               |
| Pal15        |  Myeong2017                |   Ravi      |  Shipp et al. 2018               |
| Chenab       |  Shipp et al. 2018         |   Elqui     |  Shipp et al. 2018               |
| Aliqa Uma    |  Shipp et al. 2018         |   Turbio    |  Shipp et al. 2018               |
| Willka Yaku  |  Shipp et al. 2018         | Turranburra |  Shipp et al. 2018               |
| Wambelong    |  Shipp et al. 2018         | Palca       |  Shipp et al. 2018               |
| Corvus       |  Mateu et al. 2018         | Tucana III  |  Shipp et al. 2018               |
| 20.0-1       |  Mateu et al. 2018         | Gaia-[1,5]  |  Malhan et al. 2018              |
| Phlegethon   |  Ibata et al. 2019         |             |                                  |

----------

# Module: galstreams

## MWstreams Class

The MWstreams class returns a dictionary containing a [Footprint object](#footprint-class) for each of the streams known in the MW, indexed by the stream’s name.

The class can be easily instantiated as follows:

	mwsts=galstreams.MWStreams()

This will read the stream definitions stored in the lib directory to instantiate a Footprint object appropriately for each stream. In this example, you can access the Footprint object's RA and DEC for the Orphan stream as:

  print mwsts['Orphan'].ra, mwsts['Orphan'].dec


See the Footprint Class description below for details on object attributes and methods.

## Footprint Class

This class handles a stream’s footprint as a collection of points. A Footprint object can be instantiated by passing it a name string and a pair of longitude-latitude vectors which will be interpreted as RA/DEC or l/b if the coordinate system is indicated as equatorial or galactic respectively (via de cootype keyword), e.g., as follows::

	footprint= galstreams.Footprint(ra,dec,’Amethyst’,cootype=‘equ’)

The heliocentric distance, proper motions and radial velocity can be passed at initialization as optional arguments.

The coordinates in other reference systems of interest are computed by default and set as object attributes.

As an instance of the Footprint class, a footprint object has the following default attributes:

- footprint.ra, footprint.dec    (equatorial coords)
- footprint.l, footprint.b       (galactic coords)
- footprint.cra, .cdec, .cl, .cb (footprint’s geometric center coordinates)

If Rhel, the heliocentric distance, is given the following attributes are also created:

- footprint.Rhel
- footprint.phi, .theta       (galactocentric coords, phi=0 towards the Sun)
- footprint.Rgal              (galactocentric distance)
- footprint.xhel,.yhel,.zhel  (cartesian heliocentric coords)
- footprint.x,.y,.z           (cartesian galactocentric coords)

If proper motions are given:

- footprint.pmra, .pmdec, .pmrastar  (pmrastar=pmra*cos(dec))
- footprint.pml, .pmb, .pmlstar       (pmlstar=pml*cos(dec))

If radial velocity is given:

- footprint.vrad

If all above are given:

- footprint.vxhel,.vyhel,.vzhel  (cartesian heliocentric vels)
- footprint.vx,.vy,.vz           (cartesian galactocentric vels)

A convenience method is provided to apply a mask to all array attributes of a Footprint object:

	Footprint.mask_footprint(mask)

For full details, please see the doc-string for the Footprint class.

----------

**FILES PROVIDED**

- Libraries:
	* galstreams.py
	* gcutils.py
	* bovy_coords.py
	* pyutils.py
	* lib

- Executable programs:
  * work in progress - stand-alone script to make a quick plot of the MW library in user-selected coords

- Documentation
   * README.md
