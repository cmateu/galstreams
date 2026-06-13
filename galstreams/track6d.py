import numpy as np
import astropy.table
import astropy.coordinates as ac
import astropy.units as u
import gala
import gala.coordinates as gc
import gala.dynamics as gd

class Track6D:

    def __init__(self, track_name, track_file, summary_file,
                 stream_name=None, stream_shortname=None,
                 track_reference=' ', track_discovery_references=' ',
                 verbose=True):

        """
        Track6D: A Stellar Stream's Track realization in 6D. See the list of attributes below.

        Parameters
        ==========

        - `track_name` : str

            Unique identifier for the stream's track realization. Not necesarily identical to stream_name, e.g. if
            more than one track for the same stream is available

        - `track_file` : str

            Input ecsv file containing knots to initialize 6D track stream realization

        - `summary_file` : str

            Input ecsv file containing end point, mid point and pole coordinates (6D, 6D and 3D)

        Attributes:
        ===========

        - `track_name` : str

            Unique identifier for the stream's track realization

        - `stream_name`: str

            Stream name

        - `stream_shortname`: str

            Stream short name

        - `stream_frame`: astropy coordinate frame

        - `track` : astropy.coordinates.SkyCoord Object

            Contains the track 6D info. By default initialized in icrs frame

        - `length`: astropy.Quantity Object contains the angular length measured along the track

        - `InfoFlags`: string - 4-bits indicate available (or assumed) data.

            + bit 0: 0 = great circle by construction
            + bit 1: 0 = no distance track available (only mean or central value reported)
            + bit 2: 0 = no proper motion data available (only mean or central value reported)
            + bit 3: 0 = no radial velocity data available (only mean or central value reported)

        - `end_points`: 2-element astropy.coordinates.SkyCoord Object with end point coordinates

        - `mid_point`: astropy.coordinates.SkyCoord Object with stream's mid-point coordinates (phi1=0)

        - `mid_pole`: astropy.coordinates.SkyCoord Object heliocentric pole at mid_point

        - `poly_sc`: astropy.coordinates.SkyCoord Object containing vertices for stream's polygon footprint

        - `mid_pole_gsr`: astropy.coordinates.SkyCoord Object. GSR pole at phi1=0

        - `pole_track_helio`: astropy.coordinates.SkyCoord Object heliocentric pole track (galactic coordinates by default)

        - `pole_track_gsr`: astropy.coordinates.SkyCoord Object GSR pole track (galactic coordinates by default)

        - `angular_momentum_helio`: list object with spherical components (modulus, lat, lon)
            for the angular momentum of each point along the track, computed in
            a heliocentric frame at rest w.r.t. the GSR

        WARNING: angular momentum and pole tracks have length track.size-1
        """

        #Initialize a (new) Footprint6D object

        #First the track name (track's identifier, for popular streams there can be more than one track for a given stream)
        self.track_name = str(track_name)

        #Stream's name
        if stream_name is not None:
            self.stream_name = str(stream_name)
        else:
            self.stream_name = self.track_name.copy()

        #Stream's short name
        self.stream_sname = str(stream_shortname)
        if stream_shortname is not None:
            self.stream_shortname = str(stream_shortname)
        else:
            self.stream_shortname = self.stream_name[:5]

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

        #Width attributes
        self.track_width = dict()
        for k in ['width_phi2','width_pm_phi1_cosphi2','width_pm_phi2']:
            self.track_width[k] = sfile[k][0]


        #Set up stream's coordinate frame
        if np.float64(gala.__version__[:3])<=1.4:
            self.stream_frame = gc.GreatCircleICRSFrame(pole=self.mid_pole, ra0=self.mid_point.icrs.ra)
        else:
            self.stream_frame = gc.GreatCircleICRSFrame.from_pole_ra0(
               pole=self.mid_pole,
               ra0=self.mid_point.icrs.ra,
               origin_disambiguate=self.mid_point.icrs
            )

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
        """
        Compute angular momentum for each point in the track.

        By default it returns the spherical components of the angular
        momentum in the heliocentric and galactocentric reference
        frames at rest with respect to the GSR. If `return_cartesian = True`
        it will return Cartesian components.
        """

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

        if return_cartesian:
            return L
        else:
            #Force lat,lon to be in deg because leaving them in rad is asking for trouble
            return (L[0], L[1].to(u.deg), L[2].to(u.deg) )


    def get_pole_tracks(self, use_gsr_default=True):
        """
        Compute pole at each point in the track. This is obtained by computing,
        at each point, the normal or cross product between said point and the
        contiguous point in the track
        """

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
        #m = lat<0.*u.deg
        #lon[m] = lon[m] + 180.*u.deg
        #lat[m] = np.abs(lat[m])
        pole_track_gsr = ac.SkyCoord(lon=lon, lat=lat, frame=ac.Galactocentric(), representation_type='spherical')

        return pole_track_helio, pole_track_gsr

    def create_sky_polygon_footprint_from_track(self, width=1.*u.deg, phi2_offset=0.*u.deg):

        poly_sc = create_sky_polygon_footprint_from_track(self.track, frame=self.stream_frame, width=width, phi2_offset=phi2_offset)

        return poly_sc

    def get_adql_query_from_polygon(self):

        return get_adql_query_from_polygon(self.poly_sc)

    def get_mask_in_poly_footprint(self,coo):
        """
        Return a mask for  points in input SkyCoords object that are inside polygon footprint.

        Parameters
        ==========

        - `coo` : astropy.coordinates.SkyCoord object

        Returns
        =======

        - `mask` : boolean mask array, same number of elements as coo
        """

        return get_mask_in_poly_footprint(poly_sc=self.poly_sc, coo=coo, stream_frame=self.stream_frame)



    def resample_stream_track(self, dphi1=0.02*u.deg):

        """ In construction... """


def create_sky_polygon_footprint_from_track(SkyCoordTrack, frame, width=1.*u.deg, phi2_offset=0.*u.deg):

    """
    Create the Polygon Footprint from the celestial track.
    The polygon is created by shifting the track in phi2 by a given width.

    Inputs:
    =======

    - `SkyCoordTrack`: track SkyCoord object from a MWStreams library stream (`mws[st].track`)

    - `frame`: None. Astropy coordinate frame to set up the polygon by offsetting the track by a given width.
           The default is to use the Track6D's own stream frame Track6D.stream_frame

    Parameters
    ==========

    - `phi2_offset`: astropy.Quantity object
      The offset in phi2 that will be applied to the track to create the polygon footprint (default 0)

    - `width`: astropy.Quantity object
      The total width of the polygon footprint to be created around track+phi2_offset
    """

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


def compute_angular_momentum_track(track, return_cartesian = False):

    """
    Compute angular momentum for each point in the track.
    By default it returns the spherical components of the angular momentum in
    the heliocentric and galactocentric reference frames at rest with respect to
    the GSR. Set `return_cartesian = True` to get cartesian components


    Parameters:
    =======

    - `track` : SkyCoord object

    - `return_cartesian` : If True returns cartesian coordinates. If False,
                           returns spherical coords (astropy format mod, lat, lon)

    Returns:
    ========

    - `L` : list object with compoments of angular momentum vector. By default
            returns spherical components modulus, lat, lon
    """

    tr = track.cartesian

    pos = ac.CartesianRepresentation(x = tr.x, y = tr.y, z = tr.z)
    vel = ac.CartesianDifferential(d_x = tr.differentials['s'].d_x, d_y = tr.differentials['s'].d_y, d_z = tr.differentials['s'].d_z)
    psp = gd.PhaseSpacePosition(pos=pos, vel=vel)
    L = psp.angular_momentum()

    if return_cartesian:
        return L
    else:
        L_sph = ac.cartesian_to_spherical(x = L[0], y = L[1], z = L[2])
        return L_sph
