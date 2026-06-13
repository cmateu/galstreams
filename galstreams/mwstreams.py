import os.path
import numpy as np
import astropy.table
import astropy.coordinates as ac
import astropy.units as u

from .track6d import Track6D

SUMMARY_FILE_TEMPLATE = "{tdir}/track.{imp}.{stname}.{tref}.summary.ecsv"
TRACK_FILE_TEMPLATE   = "{tdir}/track.{imp}.{stname}.{tref}.ecsv"

#---------MW Streams class--------------------------------------------------------------------------------
class MWStreams(dict):

    def __init__(self, verbose=False, implement_Off=False, print_topcat_friendly_files=False):

        #A MWStreams object is a dictionary in which each entry is a Footprint object, indexed by each stream's name.
        #There's also a mandatory summary entry, which is a Pandas DataFrame with summary attributes for the full library

        #Initialize empty dictionary
        dict.__init__(self)

        #Read in the master logs
        tdir = os.path.dirname(os.path.realpath(__file__))
        #master logs
        filepath = "{path}/{filen}".format(path=tdir+"/lib/",filen='master_log.txt')
        lmaster = astropy.table.Table.read(filepath,format='ascii.commented_header').to_pandas()
        filepath = "{path}/{filen}".format(path=tdir+"/lib/",filen='master_log.discovery_refs.txt')
        lmaster_discovery = astropy.table.Table.read(filepath,format='ascii.commented_header').to_pandas(index='Name')
        filepath = "{path}/{filen}".format(path=tdir+"/lib/",filen='master_log.comments.txt')
        lmaster_comments = astropy.table.Table.read(filepath,format='ascii.commented_header').to_pandas(index='Name')

        lmaster["On"] = lmaster["On"].astype('bool') #this attribute controls whether a given track is included or not

        #SkyCoords objects will be created for each of these dicts after the full library has been created
        attributes = ['ra','dec','distance','pm_ra_cosdec','pm_dec','radial_velocity']
        units = [u.deg, u.deg, u.kpc, u.mas/u.yr, u.mas/u.yr, u.km/u.s]
        end_o_dic = {k: np.array([])*uu  for k,uu in zip(attributes,units) }
        end_f_dic = {k: np.array([])*uu  for k,uu in zip(attributes,units) }
        mid_point_dic = {k: np.array([])*uu  for k,uu in zip(attributes,units) }
        mid_pole_dic  = {k: np.array([])*uu  for k,uu in zip(attributes,units) }
        info_flags = []
        #separate info_flags (easier to filter)
        has_empirical_track = np.array([],dtype=np.int32)
        has_D = np.array([],dtype=np.int32)
        has_pm = np.array([],dtype=np.int32)
        has_vrad = np.array([],dtype=np.int32)
        discovery_refs = []
        lengths = np.array([])#*u.deg
        track_widths = dict(phi2=np.array([]), pm_phi1_cosphi2=np.array([]), pm_phi2=np.array([]))

        print("Initializing galstreams library from master_log... ")
        nid = 1
        for ii in np.arange(lmaster.TrackRefs.size):

            #Create the names of the files containing the knots and summary attributes to initialize each stream
            format_params = dict(tdir=tdir+'/tracks', imp=lmaster.Imp[ii],
                                 stname=lmaster.Name[ii],
                                 tref=lmaster.TrackRefs[ii])
            summary_file = SUMMARY_FILE_TEMPLATE.format(**format_params)
            track_file = TRACK_FILE_TEMPLATE.format(**format_params)

            if verbose:
                print(f"Initializing Track6D {lmaster.TrackName[ii]} for {lmaster.Name[ii]}...")

            #Do the magic. The track is read and all attributes stored in the summary for all registered stream tracks.
            #Only the ones "On" are "realized" unless implement_Off == True
            track = Track6D(track_name=lmaster.TrackName[ii], stream_name=lmaster.Name[ii], track_reference=lmaster.TrackRefs[ii],
                            track_file=track_file, track_discovery_references=lmaster_discovery.loc[lmaster.Name[ii],'DiscoveryRefs'] ,
                            summary_file=summary_file)

            if implement_Off:
                self[lmaster.TrackName[ii]] = track
                self[lmaster.TrackName[ii]].ID = nid
                nid = nid+1
            else:
                if lmaster.On[ii]:
                    self[lmaster.TrackName[ii]] = track
                    self[lmaster.TrackName[ii]].ID = nid
                    nid = nid+1
                elif verbose: print(f"Skipping Off track {lmaster.TrackName[ii]}...")

            #Store summary attributes
            for k in attributes:
                end_o_dic[k] = np.append(end_o_dic[k], getattr(track.end_points, k)[0] )
            for k in attributes:
                end_f_dic[k] = np.append(end_f_dic[k], getattr(track.end_points, k)[1] )
            for k in attributes:
                mid_point_dic[k] = np.append(mid_point_dic[k], getattr(track.mid_point, k) )
            for k in attributes[:2]:
                mid_pole_dic[k]  = np.append(mid_pole_dic[k] , getattr(track.mid_pole, k) )
            for k in track_widths.keys():
                track_widths[k] = np.append(track_widths[k], track.track_width['width_'+k].value)


            info_flags.append(track.InfoFlags)
            has_empirical_track = np.append(has_empirical_track, np.int32(track.InfoFlags[0]))
            has_D               = np.append(has_D    , np.int32(track.InfoFlags[1]))
            has_pm              = np.append(has_pm   , np.int32(track.InfoFlags[2]))
            has_vrad            = np.append(has_vrad , np.int32(track.InfoFlags[3]))
            lengths = np.append(lengths, track.length.deg)
            discovery_refs = np.append(discovery_refs, lmaster_discovery.loc[lmaster.Name[ii],'DiscoveryRefs'] )


        #Add skycoord summary attributes to the library and selected cols to the summary table
        self.end_o = ac.SkyCoord(**end_o_dic)
        self.end_f = ac.SkyCoord(**end_f_dic)
        self.mid_point = ac.SkyCoord(**mid_point_dic)
        self.mid_pole  = ac.SkyCoord(ra=mid_pole_dic["ra"], dec=mid_pole_dic["dec"], frame='icrs')

        #Store master table as an attribute (inherits structure of lmaster dataframe)
        self.summary = lmaster.copy()

        #Stream Length
        self.summary["length"] = np.array(lengths)
        #End points
        self.summary["ra_o"] = end_o_dic["ra"].deg
        self.summary["dec_o"] = end_o_dic["dec"].deg
        self.summary["distance_o"] = end_o_dic["distance"].value
        self.summary["ra_f"] = end_f_dic["ra"].deg
        self.summary["dec_f"] = end_f_dic["dec"].deg
        self.summary["distance_f"] = end_f_dic["distance"].value
        #Mid point
        self.summary["ra_mid"] = mid_point_dic["ra"].deg
        self.summary["dec_mid"] = mid_point_dic["dec"].deg
        self.summary["distance_mid"] = mid_point_dic["distance"].value
        #Pole
        self.summary["ra_pole"] = mid_pole_dic["ra"].deg
        self.summary["dec_pole"] = mid_pole_dic["dec"].deg
        #Widths
        #Track widths in phi2,pm_phi1/phi2
        for k in track_widths.keys():
            mask = self.summary["width_"+k] == -1
            self.summary.loc[mask,"width_"+k] = np.round(track_widths[k][mask],decimals=2) #some streams have widths ~0.05deg

        #Info (InfoFlags and has_* columns is the same, but to have it on separate columns is more practical for filtering)
        self.summary["InfoFlags"] = np.array(info_flags)
        self.summary["has_empirical_track"] = has_empirical_track
        self.summary["has_D"]    = has_D
        self.summary["has_pm"]   = has_pm
        self.summary["has_vrad"] = has_vrad
        self.summary["DiscoveryRefs"] = discovery_refs

        #Index by TrackName
        self.summary.index=self.summary.TrackName

        #Create a numeric ID for each track
        self.summary["ID"] = ''
        for ii in self.summary.index:
            if self.summary.loc[ii,'On']:
                self.summary.loc[ii,'ID'] = self[ii].ID

        #If chosen by the user, when the library is instantiated, save in default location TOPCAT-friendly csv files with
        # the library's tracks, end-points, mid-points and summary table
        if print_topcat_friendly_files:
            self.print_topcat_friendly_compilation(output_root=f'{tdir}/tracks/galtreams.unique_streams')


    def all_unique_stream_names(self):
        """
        Returns all unique instances of the StreamNames in the library (a stream can have multiple tracks)

        Returns
        =======

        array
        """
        return np.unique(self.summary.Name[self.summary.On])

    def all_track_names(self, On_only=False):
        """
        Returns TrackNames available in the library (when `On_only=False`,
        equivalent to `MWStreams.summary['TrackName']`)

        Parameters:
        ===========

        - `On_only`: True/False
                 If True it returns only the names for the active tracks

        Returns
        =======

        array
        """

        if On_only: return self.keys()
        else:
            return np.array(self.summary.index)

    def get_track_names_for_stream(self, StreamName, On_only=False):
        """
        Find all the TrackNames for which the StreamName matches the input string (all or part of it)

        Parameters
        ==========

        - `StreamName` : str
                     Name of the stream for which to search TrackNames (or part of it)
        - `On` : book
             If True, returns only the active track name(s)

        Returns
        =======

        array : contains all the TrackNames for which the stream's name matches the input string
        """

        all_track_names = []

        if On_only:
            track_names = self.keys()
        else:
            track_names = self.summary.index

        for tn in track_names:
            if StreamName in self.summary.loc[tn,'Name']:
                all_track_names.append(tn)

        return all_track_names

    def print_topcat_friendly_compilation(self, output_root='galtreams.unique_streams'):
        import pandas as pd
        modes = ['track','end_point','mid_point']

        for mode in modes:

            full = dict(ra=np.array([]), dec=np.array([]), distance=np.array([]), pm_ra_cosdec=np.array([]), pm_dec=np.array([]),
                    radial_velocity=np.array([]), ID=np.array([]), StreamName=np.array([]), TrackName=np.array([]))

            for st in np.sort(list(self.keys())):
                if self.summary.loc[st,"On"]:

                    if 'track' in mode: x = self[st].track
                    elif 'end' in mode: x = self[st].end_points
                    elif 'mid' in mode: x = self[st].mid_point
                    N=np.size(x.ra)

                    full['ra']  = np.append(full['ra'], x.ra.deg)
                    full['dec'] = np.append(full['dec'], x.dec.deg)
                    full['distance'] = np.append(full['distance'], x.distance.value)
                    full['pm_ra_cosdec'] = np.append(full['pm_ra_cosdec'], x.pm_ra_cosdec.value)
                    full['pm_dec'] = np.append(full['pm_dec'], x.pm_dec.value)
                    full['radial_velocity'] = np.append(full['radial_velocity'], x.radial_velocity.value)
                    full['ID'] = np.append(full['ID'], self.summary.loc[st,"ID"] + np.zeros(N, dtype=int))
                    full['TrackName'] = np.append(full['TrackName'], [self[st].track_name,]*N)
                    full['StreamName']= np.append(full['StreamName'], [self[st].stream_name,]*N)

            full_pd = pd.DataFrame.from_dict(full)
            print(f"Creating TOPCAT-friendly library files:\n {output_root}.tracks/end_points/mid_points.csv")
            full_pd.to_csv(f'{output_root}.{mode}s.csv')
        #Print summary table
        print(f" {output_root}.summary.csv")
        self.summary.to_csv(f'{output_root}.summary.csv')

    def get_track_names_in_sky_window(self, lon_range, lat_range, frame, On_only=True, wrap_angle=0.*u.deg):
        """
        Get a list of track names for streams in a sky window with limits given by
        lon_range,lat_range in a given coordinate frame

        Parameters
        ==========

        - `lon_range` : np.array or list

            2-element array containing limits of sky window in "longitude" coordinate (e.g ra, l)

        - `lat_range` : np.array or list

            2-element array containing limits of sky window in "latitude" coordinate (e.g dec, b)

        - `frame` : AstroPy coordinate frame

            Coordinate frame corresponding to lon/lat_range coordinates provided above
        """

        #This is just so I can get the representation_component_names (don't know how to do it
        #without creating a frame instance, so, there, let's move on
        coo = ac.SkyCoord(np.array(lon_range), np.array(lat_range),frame=frame, unit=u.deg)
        n = dict((v,k) for k,v in coo.frame.representation_component_names.items())

        if np.any(lon_range<0.):
            wrap_angle = 180.*u.deg
            #print(f"Warning: negative longitudes - setting wrap_angle to 180. assuming this is not a mistake")

        llon = ac.Angle(lon_range).wrap_at(wrap_angle).deg
        llat = ac.Angle(lat_range).deg

        track_names = []
        for st in self.summary.TrackName:
            if On_only:
                if ~self.summary["On"][st]: continue
            #same for current track
            lon = getattr(self[st].track.transform_to(frame), n['lon']).wrap_at(wrap_angle).deg
            lat = getattr(self[st].track.transform_to(frame), n['lat']).deg

            mask =  (np.min(llon)<= lon) & (lon <= np.max(llon)) & (np.min(llat)<= lat) & (lat <= np.max(llat))

            if (mask.sum()>0): track_names.append(st)

        return track_names

    def plot_stream_compilation(self, ax=None, frame=ac.ICRS, C_attribute=None, plot_names='ID',
                                plot_colorbar = None, invert_axis=True, show_legend=True,
                                basemap_kwds = dict(projection='moll',lon_0=180., resolution='l'),
                                mlabels_kwds = dict(meridians=np.arange(0.,360.,30.), color=(0.65,0.65,0.65),linewidth=1., laxmax=90.),
                                plabels_kwds = dict(circles=np.arange(-75,75,15.), color=(0.65,0.65,0.65),linewidth=1.,
                                                    labels=[0,1,1,0], labelstyle='+/-' ),
                                scat_kwds = None,
                                annot_kwds = dict(xytext=(15,15),
                                                  textcoords='offset points',
                                                  arrowprops=dict(arrowstyle="-",color='k'),
                                                  horizontalalignment='center', verticalalignment='center'),
                                legend_kwds = dict(ncol=8,loc='center', columnspacing=0.5, handletextpad=0.1,
                                                   bbox_to_anchor=(0.5,-0.35), markerscale=3, fontsize='medium'),
                                cb_kwds = None,
                                exclude_streams=[], include_only=[], plot_On_only=False,
                                return_basemap_m = False,
                                verbose=False):
        """
        Plot a Mollweide sky projection map of the current MWStreams library
        object in the selected coordinate frame.
        Note: requires Basemap library

        Parameters
        ==========

        - `track` : SkyCoord object

        - `ax=None`

        - `frame` : Astropy astropy.coordinates.baseframe instance

            Coordinate frame to be used in sky plot

        - `C_attribute` : name of SkyCoord object attribute (in selected reference frame)
                      to pass `plt.scatter` as auxiliary column c
                      e.g. `'distance'`, `'pm_b'` if `frame=ac.Galactic`

        - `plot_names` : str ['ID','track_name','stream_name','stream_shortname']

        - `plot_colorbar`: Bool

            If C_attribute is passed, plot_colorbar=True by default

        - `invert_axis` : Bool

            Invert longitude axis, set to True by default to follow usual plotting convention for l/ra

        - `show_legend`: Bool

            Show legend at the bottom of the plot. Legend attributes can be passed via the legend_kwds dict

        - `basemap_kwds` : dict

            Keywords to instantiate Basemap projection. Default, Molweide projection

        - `mlabels_kwds`: dict  - default: `dict(meridians=np.arange(0.,360.,30.), color=(0.65,0.65,0.65),linewidth=1., laxmax=90.)`

            Meridian labelling keyword attributes to be passed to Basemap

        - `plabels_kwds`: dict - `default=dict(circles=np.arange(-75,75,15.), color=(0.65,0.65,0.65),linewidth=1.,
                         labels=[0,1,1,0], labelstyle='+/-' )`

            Parallel labelling keyword attributes to be passed to Basemap

        - `scat_kwds` : dict - default scat_kwds=dict(marker='.', s=30, alpha=0.8) [defaults change if C_attribute is passed]

            Plotting keyword attributes to be passed to plt.scatter

        - `annot_kwds` : dict

            Text and arrow attributes to be passed to annotate

        - `legend_kwds` : dict

            Legend attributes to be passed to plt.legend

        - `cb_kwds` : dict - default = `dict(label=C_attribute,  shrink=0.5)`

            Colorbar attributes to be passed to plt.colorbar

        - `exclude_streams`: list of stream TrackNames

            TrackNames for streams *not* to be included in the plot

        - `include_only`: list of stream TrackNames

            Only the TrackNames provided in this list will be plotted

        - `plot_On_only`: False

            Plot only a single (default) track for each stream (i.e. On `attribute == True`).
            The default is to plot everything in the current library object

        - `return_basemap_m` : False

            Return Basemap projection function

        - `verbose`: False

            Not doing anything right now if set to True

        Returns
        =======

        `ax` : Current axes object
        """

        if ax is None:
            fig = plt.figure(1,figsize=(17,12))
            ax = fig.add_subplot(111)


        if 'Basemap' not in sys.modules:
            from mpl_toolkits.basemap import Basemap

        m = Basemap(**basemap_kwds)

        m.drawmeridians(**mlabels_kwds)
        m.drawparallels(**plabels_kwds)
        m.drawmapboundary()

        fr = frame

        if len(include_only) > 0:
            keys_to_plot = include_only
        else:
            if plot_On_only:
                keys_to_plot = self.summary["TrackName"][self.summary.loc[:,"On"]]
                print(f"Plotting On-only streams ({len(keys_to_plot)})")
            else:
                keys_to_plot = self.keys()

        msg_flag = 0

        for st in keys_to_plot:

            if st in exclude_streams:
                continue

            #short way to ensure if plot_On_only=True plots everything, if False only On tracks are plotted

            #Get representation names for selected frame and flip dict around
            l = self[st].track.transform_to(fr).representation_component_names
            n = dict((v,k) for k,v in l.items())

            if plot_names is None:
                label, alabel = None, None
            elif 'ID' not in plot_names:
                label = "{Name}".format(Name = getattr(self[st],plot_names) )
                alabel = label
            else:
                label="{ID:.0f}={Name}".format(ID=self[st].ID,Name=self[st].track_name)
                alabel="{ID:.0f}".format(ID=self[st].ID)


            #Transform the current stream's track to selected frame
            coo = self[st].track.transform_to(fr)
            x,y = m( getattr(coo,n['lon']).value , getattr(coo, n['lat']).value )

            if scat_kwds is None:
                if C_attribute is None:
                    scat_kwds=dict(marker='.', s=30, alpha=0.8)
                elif C_attribute == 'distance':
                    scat_kwds=dict(marker='.', s=30, alpha=0.8, vmin=0., vmax=100.) #reasonable limits for distance plot
                else:
                    scat_kwds=dict(marker='.', s=30, alpha=0.8, vmin=-10., vmax=10.) #reasonable limits to plot pms

            #Extra attribute to plot
            if C_attribute is not None:
                try: c = getattr(coo, C_attribute).value
                except AttributeError:
                    c = None
                    if msg_flag==0:
                        print('WARNING: Invalid attribute selected. If not a spelling error, '
                              'you are probably trying to plot an attribute in a different '
                              'coordinate frame as the one selected. '
                              'This is currently not supported. '
                              'Plotting without C_attribute aux column for now...')
                    msg_flag = 1
            else:
                c = None

            im = ax.scatter(x,y, c=c,  **scat_kwds, label=label)


            #Using end_point to place labels
            coo = self[st].mid_point.transform_to(fr)
            xl,yl = m( getattr(coo,n['lon']).value , getattr(coo, n['lat']).value )
            xy_stream = xl, yl
            ax.annotate(alabel, xy=xy_stream,  xycoords='data', **annot_kwds)

        ax.grid(ls=':')

        if show_legend:
            ax.legend(**legend_kwds)

        if cb_kwds is None and C_attribute is not None:
            cb_kwds = dict(label=C_attribute,  shrink=0.5)

        if C_attribute is not None and plot_colorbar is None: plot_colorbar=True
        if plot_colorbar:
            plt.colorbar(im, ax=ax, **cb_kwds)

        #Follow the usual convention for Galactic and ICRS to invert the l/ra axis
        if invert_axis:
            ax.invert_xaxis()

        if return_basemap_m:
            return ax, m
        else:
            return ax
