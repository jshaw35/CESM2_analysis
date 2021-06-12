from imports import *
from functions import *

np.seterr(divide='ignore', invalid='ignore') # this maybe isn't working??

mpl.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = r"\usepackage{bm} \usepackage{amsmath}" # could I add siunitx?

class Cloud_Metric(object):
    '''
    Class for comparing model runs with GOCCP and CERES-EBAF satellite data.
    Produces latitudinal averages, comparison plots, and model bias plots 
    seasonally, and with different mapping projections.
    
    Initialize with path to directory where model output is stored.
    Cases are added with metric.add_case(casename)
    Observational data is loaded automatically.
    
    Jonah Shaw, 20200409
    '''
    
    def __init__(self, casedir):
        np.seterr(divide='ignore', invalid='ignore') # Fails to remove error message.
        for cls in reversed(self.__class__.mro()):
            if hasattr(cls, 'init'):
                cls.init(self, casedir)

    def init(self, casedir):
        # Catch directory issues from the get-go
        if not os.path.exists(casedir):
            print('ERROR: The case directory %s does not exist' % casedir)
            return None
        
        else:
            self.case_dir = casedir
            self.cases = {}
            self.case_labels = []
            self.colors = sns.color_palette("colorblind")
#             self.colors = ['red','orange','yellow','green','blue','purple','pink']
            self.__addlistsanddicts()
            
            # Prep for NCL-like plotting:
            self.contour_levels = [-6,-5,-4,-3,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3,4,5,6]
            self.color_map = cmaps.ncl_default
            
#             self.__load_GOCCP_data()
            try:
                self.__load_GOCCP_data()
                
            except:
                print('ERROR: Could not load GOCCP data.')
                print("Error: ", sys.exc_info())
                return None
            
            try:
                self.__load_CALIOP_olimpia()
            except:
                print("Could not load CALIOP SLFs.")
                self.incloud_caliop_slf = None
                self.ct_caliop_slf = None
                return None

            try:
                self.__load_CERESEBAF()
            except:
                print('Failed to load CERES-EBAF data.')
                self.ceres_data = None
                
            try:
                self.__load_ISCCP()
            except:
                print('Failed to load ISCCP data.')
                self.isccp_data = None
                
            try:
                self.__load_MISR()
            except:
                print('Failed to load MISR data.')
                self.MISR_data = None
                
            try:
                self.__load_MODIS()
            except:
                print('Failed to load MODIS data.')
                self.modis_data = None
    
    
    def __addlistsanddicts(self):
        '''
        Add some common lists and dictionaries here to declutter elsewhere.
        '''
    
        self.layer_prefixes = {'CLDLOW_CAL': ['CLDLOW_CAL','CLDLOW_CAL_LIQ',
                                              'CLDLOW_CAL_ICE','CLDLOW_CAL_UN'],
                          'CLDMED_CAL': ['CLDMED_CAL','CLDMED_CAL_LIQ',
                                         'CLDMED_CAL_ICE','CLDMED_CAL_UN'], 
                          'CLDHGH_CAL': ['CLDHGH_CAL','CLDHGH_CAL_LIQ',
                                         'CLDHGH_CAL_ICE','CLDHGH_CAL_UN'],
                          'CLDTOT_CAL': ['CLDTOT_CAL','CLDTOT_CAL_LIQ',
                                         'CLDTOT_CAL_ICE'],#,'CLDTOT_CAL_UN'],
                          'CAL_UN':     ['CLDTOT_CAL_UN','CLDLOW_CAL_UN','CLDMED_CAL_UN',
                                         'CLDHGH_CAL_UN'],
                          'CAL_LIQ':    ['CLDTOT_CAL_LIQ','CLDLOW_CAL_LIQ','CLDMED_CAL_LIQ',
                                         'CLDHGH_CAL_LIQ'],
                          'CAL_ICE':    ['CLDTOT_CAL_ICE','CLDLOW_CAL_ICE','CLDMED_CAL_ICE',
                                         'CLDHGH_CAL_ICE'],
                          'CLD_CAL':     ['CLDTOT_CAL','CLDLOW_CAL','CLDMED_CAL','CLDHGH_CAL']}
        self.seasons = ['DJF','MAM','JJA','SON']
        self.projdict = {'PlateCarree':ccrs.PlateCarree(),'Arctic':ccrs.NorthPolarStereo(),
                         'Mollweide':ccrs.Mollweide(central_longitude=180), "Antarctic":ccrs.SouthPolarStereo(),
                         'Robinson':ccrs.Robinson()}
        # Dictionary with GOCCP variables as keys and COSP variables as values: (clccalipso not in COSP?)
        self.name_dict = {
                     'cllcalipso_liq':'CLDLOW_CAL_LIQ','clmcalipso_liq':'CLDMED_CAL_LIQ',
                     'clhcalipso_liq':'CLDHGH_CAL_LIQ','cltcalipso_liq':'CLDTOT_CAL_LIQ',
                     'cllcalipso_ice':'CLDLOW_CAL_ICE','clmcalipso_ice':'CLDMED_CAL_ICE',
                     'clhcalipso_ice':'CLDHGH_CAL_ICE','cltcalipso_ice':'CLDTOT_CAL_ICE',
                     'cllcalipso':'CLDLOW_CAL','clmcalipso':'CLDMED_CAL',
                     'clhcalipso':'CLDHGH_CAL','cltcalipso':'CLDTOT_CAL',
                     'cllcalipso_un':'CLDLOW_CAL_UN','clmcalipso_un':'CLDMED_CAL_UN',
                     'clhcalipso_un':'CLDHGH_CAL_UN','cltcalipso_un':'CLDTOT_CAL_UN'
                    }
        self.ceres_var_dict = {
            'toa_sw_all_mon':'FSNT','toa_lw_all_mon':'FLNT', # FSNT is wrong (both sign and magnitude)
            'toa_sw_clr_c_mon':'FSNTC','toa_lw_clr_c_mon':'FLNTC', # _c_ means that it is just from clear footprints
            'solar_mon':'SOLIN',
#             'cldtau_total_day_mon':'CLDTAU','solar_mon':'SOLIN',
        }
#             'toa_sw_clr_t_mon':'FSNTC','toa_lw_clr_t_mon':'FLNTC', # _t_ means that it is derived from all observations and is more consistent with the model variable
#             'toa_cre_sw_mon':'SWCF','toa_cre_lw_mon':'LWCF',
# CERES variables without CAM6 analogs:
    # toa_net_all_mon
    # toa_net_clr_c_mon
    # cldarea_total_daynight_mon
    # cldpress_total_daynight_mon
    # cldtemp_total_daynight_mon
    
        self.isccp_var_dict = {
            'cltisccp':'CLDTOT_ISCCP',
        }
            
        self.modis_vars_dict = {
            'Cloud_Fraction_Retrieval_Total_Mean':'CLTMODIS',
            'Cloud_Fraction_Retrieval_Low_Mean':'CLLMODIS',
            'Cloud_Fraction_Retrieval_Mid_Mean':'CLMMODIS',
            'Cloud_Fraction_Retrieval_High_Mean':'CLHMODIS',
            'Cloud_Fraction_Retrieval_Liquid_Mean':'CLWMODIS',
            'Cloud_Fraction_Retrieval_Ice_Mean':'CLIMODIS',
            'Optical_Thickness_vs_Cloud_Top_Pressure':'CLMODIS',
            'Liquid_Path_Mean':'LWPMODIS',
            'Ice_Path_Mean':'IWPMODIS',
            'Cloud_Optical_Thickness_Liquid_Mean':'TAUWMODIS',
            'Cloud_Optical_Thickness_Liquid_MeanLog10':'TAUWLOGMODIS',
            'Cloud_Optical_Thickness_Ice_MeanLog10':'TAUILOGMODIS',
            'Cloud_Optical_Thickness_Ice_Mean':'TAUIMODIS',
            'Cloud_Optical_Thickness_Total_Mean':'TAUTMODIS',
            'Cloud_Optical_Thickness_Total_MeanLog10':'TAUTLOGMODIS',
            'Cloud_Top_Pressure_Total_Mean':'PCTMODIS',
            'Cloud_Particle_Size_Liquid_Mean':'REFFCLWMODIS',
            'Cloud_Particle_Size_Ice_Mean':'REFFCLIMODIS',
            'Cloud_Optical_Thickness':'cosp_tau_modis',
            'Cloud_Top_Pressure':'cosp_prs',
        }
        
        self.modis_cloud_fracs = ['CLTMODIS','CLLMODIS','CLMMODIS',
                                  'CLHMODIS','CLWMODIS','CLIMODIS',
                                  'CLMODIS']
            
        self.lat_bounds = [[-82,-70],[-70,-60],[-60,-50],[-50,-40],[-40,-30],
              [-30,-20],[-20,-10],[-10,0],[0,10],[10,20],[20,30],
              [30,40],[40,50],[50,60],[60,70],[70,82]]
        
        self.seas_dict = {"DJF":0,"JJA":1,"MAM":2,"SON":3}
        
        self.var_label_dict = {'CLDTOT_CAL':'Total Cloud Fraction (\%)', # had a '\n'
                               'CLDTOT_CAL_LIQ':'Liquid Cloud \n Fraction (\%)',
                               'CLDTOT_CAL_ICE':'Ice Cloud \n Fraction (\%)',
                               'CLDTOT_ISCCP':'Total Cloud Fraction (\%)',
                               'CLDTOT_MISR':'Total Cloud Fraction (\%)',
                              }
        
        self.months = ['J','F','M','A','M','J','J','A','S','O','N','D'] # month initials
        
        
        
    def __load_GOCCP_data(self):
        '''
        Load GOCCP data for model comparison. Relabel variables to match model output.
        
        Original data source:
        ftp://ftp.climserv.ipsl.polytechnique.fr/cfmip/GOCCP_v3/2D_Maps/grid_2x2xL40/
        '''
        print('Loading GOCCP data...', end = '')
        
        processed_path = '/glade/u/home/jonahshaw/w/obs/CALIPSO/GOCCP/2Ddata/1.25x0.9_python_interp/amip_processed.nc'
        if os.path.exists(processed_path):
            self.goccp_data = xr.open_dataset(processed_path)
            print('done.')
            return
#         goccp_dir = '/glade/u/home/jonahshaw/w/obs/CALIPSO/GOCCP/2Ddata/1.25x0.9_interpolation/'
        goccp_dir = '/glade/u/home/jonahshaw/w/obs/CALIPSO/GOCCP/2Ddata/1.25x0.9_python_interp/'
        #'/nird/home/jonahks/p/jonahks/GOCCP_data/'
                
        # These have an extra dimenson and need to indexed at 0
        unvars = ['cllcalipso_un','clmcalipso_un','clhcalipso_un','cltcalipso_un']
        
        if not os.path.exists(goccp_dir):
            print('Could not find GOCCP directory %s' % goccp_dir)
            return
        
        # Assemble lists for each file type
        # Assume that we are using all available files (no need to select by data)
        
        phase = [] # netCDF files with cloud phase data
        cloud = []   # netCDF files with cloud amount data
        
        for year in os.listdir(goccp_dir): # iterate over years (2009-2013)
            _prepath = '%s%s/' % (goccp_dir, year)
            _files = os.listdir(_prepath)
            
            _phasefiles = [_prepath + x for x in _files if 'Phase' in x]
            _cloudfiles = [_prepath + x for x in _files if 'MapLowMidHigh330m' in x]
            
            phase = phase + _phasefiles
            cloud = cloud + _cloudfiles
            
        phase.sort()
        cloud.sort()
        
        _phase_ds = xr.open_mfdataset(phase, combine='by_coords')
        _cloud_ds = xr.open_mfdataset(cloud, combine='by_coords')
        
        _goccp_data = xr.merge([_phase_ds, _cloud_ds])

        # Select right value for undefined clouds.
        for i in unvars:
             _goccp_data[i] = _goccp_data[i].isel(cat1=0)
        # Scale cloud fraction values to a percent value to match COSP
        goccp_vars = list(self.name_dict.keys())

        for i in goccp_vars: #think maybe: with xr.set_options(keep_attrs=True)
            _goccp_data[i] = 100*_goccp_data[i]

        # Quickly add a variable to generate weights
        _goccp_data = add_weights(_goccp_data) # this is causing an issue before I interpolate. Error:  (<class 'KeyError'>, KeyError('lat'), <traceback object at 0x2b3808c7beb0>)

        # Rename variables so they will match COSP names and index correctly
        self.goccp_data = _goccp_data.rename(self.name_dict)
                
        print('done.')
        
    def __load_CALIOP_olimpia(self):
        print('Loading CALIOP SLFs...', end = '')
        try:
            self.ct_caliop_slf = xr.open_dataset('caliop_olimpia/ct_slf_olimpia/cloudtop_slfs.nc')
        except:
            print('Could not load cloudtop CALIOP slfs from caliop_olimpia/ct_slf_olimpia/cloudtop_slfs.nc')
            return
        try:
            self.incloud_caliop_slf = xr.open_dataset('caliop_olimpia/incloud_slf_olimpia/incloud_slfs.nc')
        except: 
            print('Could not load incloud CALIOP slfs from caliop_olimpia/incloud_slf_olimpia/incloud_slfs.nc')
            return
            
        print("done")
        
    def __load_CERESEBAF(self):
        '''
        Load interpolated CERES-EBAF data and rename variables for easier comparison/incorporation.
        '''
        print('Loading CERES-EBAF fluxes...', end = '')
#         _ceres_data = xr.open_dataset('/glade/u/home/jonahshaw/obs/CERES_EBAF/CERES_EBAF-TOA_Ed4.1_Subset_200003-202002.nc')
        _ceres_data = xr.open_dataset('/glade/work/jonahshaw/obs/CERES/CERES_EBAF-TOA_Ed4.1_Subset_200003-202102.nc')
#         _ceres_data = _ceres_data.sel(time=slice('2000-01-01', '2015-12-31')) # This was a problem.
        _ceres_data = add_weights(_ceres_data)
        # Flip FLNSC so it matches model convention (net-LW is down, net-SW is up)
#         _ceres_data['sfc_net_lw_clr_t_mon'] = -1*_ceres_data['sfc_net_lw_clr_t_mon']
        
        # Derive Cloud Radiative Effect Variables
        _ceres_data['SWCF'] = (_ceres_data['toa_sw_clr_c_mon'] - _ceres_data['toa_sw_all_mon']).assign_attrs(
            {'units': 'W/m2','long_name': 'Shortwave cloud forcing'})
        _ceres_data['LWCF'] = (_ceres_data['toa_lw_clr_c_mon'] - _ceres_data['toa_lw_all_mon']).assign_attrs(
            {'units': 'W/m2','long_name': 'Longwave cloud forcing'})
        # Rename variables so they will match standard CAM variables names and index correctly
        self.ceres_data = _ceres_data.rename(self.ceres_var_dict)

        print('done.')
        
    def __load_ISCCP(self):
        '''
        Load ISCCP cloud total processed for COSP ISCCP simulator.
        '''
        print('Loading ISCCP cloud total...', end = '')
        _isccp_data = xr.open_dataset('/glade/u/home/jonahshaw/w/obs/ISCCP/cltisccp_198307-200806.nc')
        _isccp_data = _isccp_data.rename(self.isccp_var_dict)
        
        self.isccp_data = _isccp_data
        
        print('done.')
        
    def __load_MISR(self):
        '''
        Load monthly averaged MISR cldtau-cldheight histograms from 2000/04-2020/05.
        '''
        print('Loading MISR cloud histograms...', end = '')
        
        misr_dir = '/glade/work/jonahshaw/obs/MISR/interp_0.9x1.25/'
        misr_monthly_file = 'clMISR_obs4MIPs_MISR_V7_monthly_averages.nc'
        misr_monthly = xr.open_dataset('%s/%s' % (misr_dir,misr_monthly_file))
        
        misr_monthly = misr_monthly.assign_coords({'cth':misr_monthly.cth*1e-3}) # convert cloudtop height to km to match COSP
        misr_monthly['cth'] = misr_monthly['cth'].assign_attrs({'units' : 'km'}) # rename variable to km

        misr_monthly = misr_monthly.rename({'clMISR':'CLD_MISR','tau':'cosp_tau','cth':'cosp_htmisr'})

        # Assign fake time attributes so it will index correctly with the metric
        misr_monthly = misr_monthly.rename({'month':'time'})
        misr_monthly['time'].attrs['calendar'] = '360_day'
        misr_monthly['time'].attrs['units'] = "months since 1999-12-01"

        misr_monthly = xr.decode_cf(misr_monthly)
        
        self.misr_data = misr_monthly
        
        print('done.')
        
    def __load_MODIS(self):
        '''
        Load MODIS cldtau-cldheight histogram climatology from 2002/07 to 2010/07.
        '''
        
        print('Loading MODIS cloud histograms...', end = '')
        
        modis_dir = '/glade/work/jonahshaw/obs/MODIS/'
        modis_climo_file = 'MCD03_M3_NC_200207to201007.V01.nc'
        modis_climo = xr.open_dataset('%s/%s' % (modis_dir,modis_climo_file))
        
        # convert cloudtop pressure from Pa to mb to match COSP
        modis_climo = modis_climo.assign_coords({'Cloud_Top_Pressure':modis_climo.Cloud_Top_Pressure*1e-2}) 
        modis_climo['Cloud_Top_Pressure'] = modis_climo['Cloud_Top_Pressure'].assign_attrs({'units' : 'mb'}) # rename variable to km
        
        modis_climo = modis_climo.rename(self.modis_vars_dict)
        
        for i in self.modis_cloud_fracs: # Convert from cloud fraction to cloud percent
            modis_climo[i] = 100*modis_climo[i]
        
        self.modis_data = modis_climo
        
        print('done.')
        
        
    def add_case(self, case, path=None, label=None, t_range=None): # Dictionary can be overwritten
        # Add a Model_case object to the cases dictionary
        if path == None: # Default argument values are evaluated at function define-time, but self is an argument only available at function call time.
            path = self.case_dir
        self.cases[case] = Model_case(path, case, t_range=t_range) # opportunity for kwargs here
        if label:
            self.cases[case].label = label
            self.case_labels += [label]
        else:
            self.case_labels += [self.cases[case].label]
    
    def get_case(self, case):
        if case not in self.cases: # check if it is already in the dictionary
            print('Could not find %s in the SLF_Metric object' % case)
            return None
        else:
            return self.cases[case] # should have a check for the object type
    
    def get_cases(self):
        return self.cases
    
    def plot1D(self, var, layers=False, seasonal=False,season=None, bias=False, lat_range=None):
        '''
        General line plot function averaging over longitudes.
        '''
        # Handle input errors
        if (seasonal and season):
            print('Cannot plot both season %s and seasonal: True. \n')
            return None     
        if (season and season not in self.seasons):
            print('%s not a valid argument for "season"' % season, self.seasons)
            return None
        
        # call appropriate plotting wrapper
        if layers:
            try: 
                varlist = self.layer_prefixes[var]
                self.__check_for_vars(varlist)
            except:
                print('''Layers prefix %s not found. Select from: \n
                %s ''' % (var, str(list(self.layer_prefixes.keys()))))
                return None
            if seasonal:
                return self.__seasonallayers1Dplotwrapper(varlist, bias=bias,
                            lat_range=lat_range)
            else:    
                return self.__layers1Dplotwrapper(varlist, bias=bias, 
                            lat_range=lat_range, season=season)
        self.__check_for_vars(var)
        if seasonal:
            return self.__seasonal1Dplotwrapper(var, bias=bias, 
                        lat_range=lat_range)
        else:
            return self.__standard1Dplotwrapper(var, bias=bias, 
                        lat_range=lat_range, season=season)
        

    def __standard1Dplotwrapper(self, var, bias=False, lat_range=None, season=None, **kwargs):
        '''Plot variable against latitude for observations and loaded model runs.'''

        # Using plt.subplots only for consistency here.
        fig, axes = plt.subplots(nrows=1,ncols=1,figsize=[6,4])
        fig.set_dpi(200)
        obs_source, obs_label = self.__get_data_source(var)
#         obs_label = obs_label.replace('_',' ') # cleans underscores from labels so Latex doesn't bug out JKS
        obs_da = obs_source
        if not bias:
            if season:
                _season_da = season_mean(obs_source[var]).where(np.absolute(obs_source['lat'])<82)
                obs_da = _season_da[self.seas_dict[season]].to_dataset(name=var)
            self.__standard1Dplot(var, obs_da, axes, bias=False, lat_range=lat_range, 
                                  label=obs_label,color='black',**kwargs)
#             self.__standard1Dplot(var, obs_source, axes,
#                                   lat_range=lat_range, label=obs_label)
        
        for k,color in zip(self.cases,self.colors):
            _run = self.cases[k]
            _da = _run.case_da
            if season:
                _season_da = season_mean(_da[var]).where(np.absolute(_da['lat'])<82)
                _da = _season_da[self.seas_dict[season]].to_dataset(name=var)
            
            # Works with percents only
            self.__standard1Dplot(var, _da, axes, bias=bias, lat_range=lat_range, label=_run.label,color=color,**kwargs)
            
        try: 
            axes.set_ylabel('%s (%s)' % (_da[var].long_name,_da[var].units))
        except:
            axes.set_ylabel('%s' % (var))
#         fig.legend()
        
        return fig

    def __layers1Dplotwrapper(self, varlist, bias=False, lat_range=None, season=None, **kwargs):
        '''Plot variable against latitude for observations and loaded model runs.'''

#         fig, axes = plt.subplots(nrows=len(varlist),ncols=1,figsize=[10,2*len(varlist)])

        fig, axes = plt.subplots(nrows=1,ncols=len(varlist),figsize=[15,5])
        obs_source, obs_label = self.__get_data_source(varlist[0])
        
        if bias:
            labels = self.case_labels
        else:
            labels = [obs_label] + self.case_labels
                
        obs_da = obs_source
        for var,ax in zip(varlist,axes):                
            if not bias:
                if season:
                    _season_da = season_mean(obs_source[var]).where(np.absolute(obs_source['lat'])<82)
                    obs_da = _season_da[self.seas_dict[season]].to_dataset(name=var)
                self.__standard1Dplot(var, obs_da, ax, bias=False, lat_range=lat_range, 
                                      label=obs_label,color='black')
            
            for k,color in zip(self.cases,self.colors):
                _run = self.cases[k]
                _da = _run.case_da
                if season:
                    _season_da = season_mean(_da[var]).where(np.absolute(_da['lat'])<82)
                    _da = _season_da[self.seas_dict[season]].to_dataset(name=var)
                
                self.__standard1Dplot(var, _da, ax, bias=bias,
                                      lat_range=lat_range,label=_run.label,color=color)
        
        fig.legend(labels=labels)
        
        xlabels = []
        for i in varlist: # jks improve labelling with self.var_label_dict
            try:
                xlabels.append(self.var_label_dict[i])
            except:
                xlabels.append(i)
        xax = axes
        
        if bias:
            self.add_labels(xlabels=xlabels, xaxes=xax, ylabels=["Bias (Model-Observation)"], 
                            yaxes=[axes[0]])
        else:
            self.add_labels(xlabels=xlabels, xaxes=xax) # JKS testing error
        self.share_ylims(axes)
        
        fig.subplots_adjust(hspace=0.5)
        
        return fig
        
    def __seasonal1Dplotwrapper(self, var, bias=False, lat_range=None, **kwargs):
        '''Plot variable against latitude for observations and loaded model runs.'''
        obs_source, obs_label = self.__get_data_source(var)
        fig, axes = plt.subplots(nrows=1,ncols=4,figsize=[15,4])
        fig.set_dpi(200)
    
        if bias: 
            labels = self.case_labels
        else:
            labels = [obs_label] + self.case_labels
            _season_da = season_mean(obs_source[var]).where(np.absolute(obs_source['lat'])<82)
            for season,ax in zip(self.seasons,axes):
                seas = _season_da.sel(season=season)
                self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
                                      bias=False, lat_range=lat_range,color='black')
            
#             for seas,ax in zip(_season_da,axes): # jks
#                 self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
#                                       bias=False, lat_range=lat_range)
            
        # Iterate cases first, and then seasons within that
        for k,color in zip(self.cases,self.colors):
            _run = self.cases[k]
            _da = _run.case_da
            
            _season_da = season_mean(_da[var])
            for season,ax in zip(self.seasons,axes):
                seas = _season_da.sel(season=season)
#                 print(seas.season.values)
                self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
                                      bias=bias,lat_range=lat_range, label=_run.label,color=color)
                
#             for seas,ax in zip(_season_da,axes): # jks
#                 self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
#                                       bias=bias,lat_range=lat_range, label=_run.label)
        
        fig.legend(labels=labels)
        
        ylabels = [var]; yax = [axes[0]]
        xlabels = self.seasons; xax = axes
#         xlabels = _season_da.coords['season'].values; xax = axes
        self.add_labels(xlabels=xlabels, xaxes=xax, ylabels=ylabels, yaxes=yax)
        self.share_ylims(axes)

        return fig

    def __seasonallayers1Dplotwrapper(self, varlist, bias=False, lat_range=None, **kwargs):
        '''Plot variables by season against latitude for observations and loaded model
        runs.'''

        fig, axes = plt.subplots(nrows=len(varlist),ncols=4,figsize=[15,3*len(varlist)])
        fig.set_dpi(200)
        obs_source, obs_label = self.__get_data_source(varlist[0])
        
        # Preprocessing step:
        for var, xax in zip(varlist, axes): # iterate over rows (i.e. layer vars)
            var_processed = []
            if bias: 
                labels = self.case_labels
            else:
                labels = [obs_label] + self.case_labels
                _season_da = season_mean(obs_source[var]).where(np.absolute(obs_source['lat'])<82)

                for season,ax in zip(self.seasons,xax):
                    seas = _season_da.sel(season=season)
                    self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
                                          bias=False, lat_range=lat_range, label=obs_label,color='black')
                
#                 for seas,ax in zip(_season_da,xax): #iterate over columns (seasons)
#                     self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
#                                           bias=False, lat_range=lat_range, label=obs_label)
            
            for k,color in zip(self.cases,self.colors):
                _run = self.cases[k]
                _da = _run.case_da
                _season_da = season_mean(_da[var])
                
                for season,ax in zip(self.seasons,xax):
                    seas = _season_da.sel(season=season)
                    self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
                                          bias=bias, lat_range=lat_range, label=_run.label,color=color)
                
#                 for seas,ax in zip(_season_da,xax): #iterate over columns (seasons)
#                     self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
#                                       bias=bias, lat_range=lat_range, label=_run.label)
                    
            self.share_ylims(xax) # Share ylims only across layers
        
        fig.legend(labels=labels)
        
        ylabels = varlist; yax = axes[:,0]
        xlabels = self.seasons; xax = axes[0,:] # jks
        self.add_labels(xlabels=xlabels, xaxes=xax, ylabels=ylabels, yaxes=yax)
        
        return fig
        
    def __standard1Dplot(self, var, da, ax, bias=False, lat_range=None, **kwargs):
        '''
        Simple 1D plotting function. Need to build in better labels+bias.'''
        obs_source, obs_label = self.__get_data_source(var)

        lat_lims = [-90,90]
        if lat_range:
            lat_lims = lat_range
        if (bias and lat_lims[0] < -82): lat_lims[0] = -82
        if (bias and lat_lims[1] > 82): lat_lims[1] = 82
            
        if 'time' not in da.dims: # catches a bug with seasonal averages
            val = da[var].mean(dim='lon', skipna=True)
        else:
            val = da[var].mean(dim = ['time','lon'], skipna=True)

        if bias:
            try: # get season observations if passed da is seasonally processed
                season = val['season'].values
                obs_period = season_mean(obs_source[var]).sel(season=season).where(np.absolute(obs_source['lat'])<82)
                obs = obs_period.mean(dim = ['lon'], skipna=True)
            except:
                obs = obs_source[var].mean(dim = ['time','lon'], skipna=True)
            
            val = val.interp_like(obs_source[var]) # quick interp fix for weird grid mismatch (bad?)
            val = (val - obs).sel(lat=slice(lat_lims[0],lat_lims[1]))

#             im = val.sel(lat=slice(lat_lims[0],lat_lims[1])).plot(ax=ax, add_legend=False, **kwargs)
            
            ax.plot(np.sin(np.pi/180*val['lat']),val,**kwargs)
            ax.set_xticks([-1,-1*np.sqrt(3)/2,-0.5,0,0.5,np.sqrt(3)/2,1])
            ax.set_xticklabels(['90S','60S','30S','0','30N','60N','90N'])
            
            ax.hlines(0, lat_lims[0], lat_lims[1], colors='gray', linestyles='dashed', label='')
        else:
#             im = val.sel(lat=slice(lat_lims[0],lat_lims[1])).plot(ax=ax, add_legend=False, **kwargs)
            val = val.sel(lat=slice(lat_lims[0],lat_lims[1]))
            
            ax.plot(np.sin(np.pi/180*val['lat']),val,**kwargs)
            ax.set_xticks([-1,-1/np.sqrt(2),-0.5,0,0.5,1/np.sqrt(2),1])
            ax.set_xticklabels(['90S','45S','30S','0','30N','45N','90N'])
#             ax.set_xticks([-1,-1*np.sqrt(3)/2,-0.5,0,0.5,np.sqrt(3)/2,1])
#             ax.set_xticklabels(['90S','60S','30S','0','30N','60N','90N'])
            
#         ax.set_xlim([-1.1,1.1]) # weird fix??
        ax.set_xlim([-1,1]) # weird fix??
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_title('')
        
    def plot2D(self, var, projection='PlateCarree',layers=False,
               seasonal=False, season=None, bias=False,**kwargs):
        '''
        General surface plot function. Calls more specific plotting functions.
        Supports projections:
        PlateCarree, NorthPolarStereo, Mollweide, SouthPolarStereo, Robinson
        ! currently avoid Mollweide because of high-latitude bug:
        https://github.com/SciTools/cartopy/issues/1604
        '''
                
        # Series of comparison to select specific plotting function to use:
        if (layers and seasonal):
            print(""""Cannot plot along 3 dimensions (case, layer, season). \n
            Please select either layers OR seasonal""")
            return None
        if (seasonal and season):
            print('Cannot plot both season %s and seasonal: True. \n')
            return None            
        if projection not in self.projdict:
            print('Did not recognize projection string.')
            print('Please select from: %s' % self.projdict.keys())
            return None            
        else: proj = self.projdict[projection]
        if (season and season not in self.seasons):
            print('%s not a valid argument for "season"' % season, self.seasons)
            return None
        
        # call appropriate plotting wrapper
        if layers:
            try: 
                varlist = self.layer_prefixes[var]
                self.__check_for_vars(varlist)
            except:
                print('''Layers prefix %s not found. Select from: \n
                %s ''' % (var, str(list(self.layer_prefixes.keys()))))
                return None
            out = self.__layersplotwrapper(varlist, proj, bias=bias, season=season,**kwargs)
        elif seasonal:
            self.__check_for_vars(var)
            out = self.__seasonalplotwrapper(var, proj, bias=bias,**kwargs)
        else:
            self.__check_for_vars(var)
            out = self.__standard2Dplotwrapper(var, proj, bias=bias, season=season,**kwargs)
                
        return out
        
    def __standard2Dplotwrapper(self, var, projection, bias=False, season=None, label=True, **kwargs):
        '''
        Create subplots object and iterate over cases to plot with correct projection, etc.
        '''
        obs_source, obs_label = self.__get_data_source(var)

        out = False
        if 'ax' in kwargs.keys(): # use user supplied axes, must be of appropriate size and just a list
            out = True
            if bias:
                case_axes = kwargs['ax']
            else:
                case_axes = kwargs['ax'][1:]
                axes = np.array(kwargs['ax'])
            del kwargs['ax'] # hmm
        
        else:
            if bias:
                fig, axes = sp_map(nrows=len(self.cases), ncols=1,
                               projection=projection, figsize=[15,2*(len(self.cases))])
                fig.set_dpi(200)
                ylabels = self.case_labels
                case_axes = axes
                if len(self.cases) == 1:  # jks handle axes if not an interable
                    case_axes = [axes]

            else:
                fig, axes = sp_map(nrows=len(self.cases)+1, ncols=1,
                               projection=projection, figsize=[20,2*(len(self.cases)+1)])
                fig.set_dpi(200)
                ylabels = [obs_label] + self.case_labels
                if len(self.cases) == 0:  # jks handle axes if not an interable
                    axes = np.array([axes])
                    case_axes = []
                else:
    #             _axes = axes if len(self.cases) == 0 else axes.flat[1:]
    #             if len(self.cases) == 0: _axes = axes
                    case_axes = axes.flat[1:]
        
#         if 'ax' in kwargs.keys(): # use user supplied axes, must be of appropriate size and just a list
#             if bias:
#                 case_axes = kwargs['ax']
#             else:
#                 case_axes = kwargs['ax'][1:]
#                 axes = np.array(kwargs['ax'])
#             del kwargs['ax']
        
        if not bias: # weird re-org here
            # jks handle season
            _da = obs_source
            _label = '%s %s' % (obs_label, self.var_label_dict[var])
            if season:
                _season_da = season_mean(obs_source[var]).where(np.absolute(obs_source['lat'])<82)
                _da = _season_da[self.seas_dict[season]].to_dataset(name=var)
            _ax,_im = self.standard2Dplot(var, _da, axes.flat[0],
                                projection=projection, bias=False, label=_label, **kwargs) # JKS include var name here
#             _axes = axes.flat[1:]
        
        for ax, k, _label in zip(case_axes, self.cases, self.case_labels):
            if bias: _label += ' Bias'
            _run = self.cases[k]
            _da = _run.case_da
            
            # Works with percents only
            if season:
                _season_da = season_mean(_da[var])
                _da = _season_da[self.seas_dict[season]].to_dataset(name=var)
                
            _ax, _im = self.standard2Dplot(var,_da, ax,projection=projection,
                                           bias=bias, label=_label, **kwargs)
                        
        if label:
            yax = axes
            try:
                xlabels = [var]; xax = [axes[0]]
            except:
                xlabels = [var]; xax = [axes]
            self.add_labels(ylabels=ylabels, yaxes=yax, xlabels=xlabels, xaxes=xax,height=1.25)
        
        if not 'contour' in kwargs.keys():
            self.share_clims(fig)
        
        if 'add_colorbar' in kwargs.keys():
            if kwargs['add_colorbar']:
                try:
                    cbar = fig.colorbar(_im, ax=axes.ravel().tolist()) # this doesn't get the right bounds
                except:
                    cbar = fig.colorbar(_im, ax=axes)
            
                if bias:
                    cbar.set_label("Bias (Model - Observations) of %s (%s)" % (_da[var].long_name,_da[var].units))
        #             cbar.set_label("Bias (%s)" % _da[var].units)
                else:
                    cbar.set_label("%s (%s)" % (_da[var].long_name,_da[var].units))#, fontsize=12)

        if out == True:
            return None,_im
        else:
            return fig, _im

    def __layersplotwrapper(self, varlist, projection, bias=False, season=None, **kwargs):
        '''
        Create subplots object and iterate over cases to plot with correct projection, etc.
        '''
        obs_source, obs_label = self.__get_data_source(varlist[0]) # layers must share a source
        
        if bias:
            fig, axes = sp_map(nrows=len(self.cases), ncols=len(varlist),
                           projection=projection, figsize=[5.25*len(varlist),3*(len(self.cases))])
            ylabels = self.case_labels
            
        else:
            fig, axes = sp_map(nrows=len(self.cases)+1, ncols=len(varlist),
                           projection=projection, figsize=[5.25*len(varlist),3*(len(self.cases)+1)])
            ylabels = [obs_label] + self.case_labels
#         fig.set_dpi(100)
            
        obs_da = obs_source
        for var,varax in zip(varlist,axes.transpose()):
            # Plot observational data, always the raw data (not bias)
            if bias:
                case_ax = varax
            else:
                if season:
                    _season_da = season_mean(obs_source[var]).where(np.absolute(obs_source['lat'])<82)
                    obs_da = _season_da[self.seas_dict[season]].to_dataset(name=var)
                    
                _ax, _im = self.standard2Dplot(var, obs_da, varax[0],
                                projection=projection, bias=False, label=obs_label, **kwargs)
                case_ax = varax[1:]
                
        
            for ax, k, _label in zip(case_ax, self.cases, self.case_labels):
                _run = self.cases[k]
                _da = _run.case_da
            
                if season:
                    _season_da = season_mean(_da[var])
                    _da = _season_da[self.seas_dict[season]].to_dataset(name=var)

                # Works with percents only
                _ax, _im = self.standard2Dplot(var,_da, ax,projection=projection,
                                               bias=bias, label=_label, **kwargs)
        
        yax = axes[:,0]
        xlabels = []
        for i in varlist: # jks improve labelling with self.var_label_dict
            try:
                xlabels.append(self.var_label_dict[i])
            except:
                xlabels.append(i)
                
        xax = axes[0,:]
        self.add_labels(ylabels=ylabels, yaxes=yax, xlabels=xlabels, xaxes=xax,height=1.25)

        if not 'contour' in kwargs.keys():
            self.share_clims(fig)
        plt.subplots_adjust(wspace=0.1, hspace=0.1)
        
        cbar = fig.colorbar(_im, ax=axes.ravel().tolist()) 
        # ^does this work after the _im object has been modified by share_clims? 
        # Yes, it does!
        if bias:
            cbar.set_label("Bias (Model - Observations)")
            
        return fig, cbar
    
    def __seasonalplotwrapper(self, var, projection, bias=False, **kwargs):
        '''
        Create subplots object and iterate over cases to plot with correct projection, etc.
        '''        
        obs_source, obs_label = self.__get_data_source(var)
        if bias:
            fig, axes = sp_map(nrows=len(self.cases), ncols=4,
                           projection=projection, figsize=[10,2*len(self.cases)])
#             fig.set_dpi(100)
            ylabels = self.case_labels
            _axes = axes
            
        else:
            fig, axes = sp_map(nrows=len(self.cases)+1, ncols=4,
                           projection=projection, figsize=[10,2*(len(self.cases)+1)])
#             fig.set_dpi(100)
            ylabels = [obs_label] + self.case_labels
            
            _season_da = season_mean(obs_source[var]).where(np.absolute(obs_source['lat'])<82)

            for season,ax in zip(self.seasons,axes[0]):
                seas = _season_da.sel(season=season)
                _ax, _im = self.standard2Dplot(var, seas.to_dataset(name=var), ax,
                                    projection=projection, bias=False,label=obs_label, **kwargs)
                
#             for seas,ax in zip(_season_da,axes[0]):
#                 _ax, _im = self.standard2Dplot(var, seas.to_dataset(name=var), ax,
#                                     projection=projection, bias=False)
            _axes = axes[1:]
            
        # Iterate cases first, and then seasons within that
        for xaxes, k, _label in zip(_axes, self.cases, self.case_labels):
            _run = self.cases[k]
            _da = _run.case_da
            
            _season_da = season_mean(_da[var])
            
            for season,ax in zip(self.seasons,xaxes):
                seas = _season_da.sel(season=season)
                _ax, _im = self.standard2Dplot(var, seas.to_dataset(name=var), ax,
                                projection=projection, bias=bias, label=_label, **kwargs)
            
#             for seas,ax in zip(_season_da,xaxes):
#                 _ax, _im = self.standard2Dplot(var, seas.to_dataset(name=var), ax,
#                                 projection=projection, bias=bias)

        if not 'contour' in kwargs.keys():
            self.share_clims(fig)
        
        cbar = fig.colorbar(_im, ax=axes.ravel().tolist())
        if bias:
            cbar.set_label("Bias (%s)" % _da[var].units)
        else:
            cbar.set_label("%s (%s)" % (_da[var].long_name,_da[var].units), fontsize=12)
        
        yax = axes[:,0]
        xlabels = self.seasons; xax = axes[0,:]
        self.add_labels(ylabels=ylabels, yaxes=yax, xlabels=xlabels, xaxes=xax,height=1.25)
                
        return fig
    
    def standard2Dplot(self, var, da, ax, projection=None, bias=False, contour=False, afunc=None,
                       label=None, **kwargs):
        '''
        Basic workhorse plotting function. needs to handle bias and seasons
        '''
        obs_source, obs_label = self.__get_data_source(var)
        
        bold_label = bfize(label)
        if len(label) > 20 or bias: # large titles need to be split
            bold_label = bold_label + ' \n '   
        
        lat_lims = [-90,90]
        if projection == ccrs.NorthPolarStereo(): 
            lat_lims = [59.5,90]
            polarCentral_set_latlim(lat_lims, ax)
        if projection == ccrs.SouthPolarStereo(): 
            lat_lims = [-90,-60]
            polarCentral_set_latlim(lat_lims, ax)
        if 'time' not in da[var].dims: # catches a bug with seasonal averages
            val = da[var].where(da['lat'] > lat_lims[0])
        else:
            val = da[var].mean(dim = 'time', skipna=True).where(
                                da['lat'] > lat_lims[0])            
                    # apply an arbitrary function (must return 2d output):
        if afunc:
            val = afunc(val)
            
        standards = {'levels':20 if contour else None, 
                     'extend':'both', 
                     'robust':True,
                     'add_colorbar':False,
                     'transform':ccrs.PlateCarree(),
                     'cmap':self.color_map if contour else 'vlag' if bias else 'viridis'} # handle 3(4) cases
#                      'cmap':'vlag' if bias else self.color_map} # handle 3(4) cases
        if ('levels' not in kwargs.keys()) or np.any(kwargs['levels'] == None):
            lvl_bool = True
        else: lvl_bool = False
        for i in standards: # adopt standards if the argument is unspecified
            if i not in kwargs.keys():
                kwargs[i] = standards[i]            
            
        if bias:
            try: # get season observations if passed da is seasonally processed
                season = val['season'].values
                obs = season_mean(obs_source[var]).sel(season=season).where(
                                  (da['lat'] > lat_lims[0]) & (np.absolute(obs_source['lat'])<82))

            except:
                if 'time' in obs_source[var].dims:
                    obs = obs_source[var].mean(dim = 'time', skipna=True).where(
                                                obs_source['lat'] > lat_lims[0]) # obs_source was da, not sure why
                else:
                    obs = obs_source[var].where(obs_source['lat'] > lat_lims[0])
                    
            # apply an arbitrary function (must return 2d output):
            if afunc:
                obs = afunc(obs)
            obs = obs.interp_like(val) # quick interp fix for weird grid mismatch (bad.)

            # Calculate global average error
            weights = add_weights(val)['cell_weight']
            val_avg = masked_average(val,weights=weights,dim=['lat','lon'])
            obs_avg = masked_average(obs,weights=weights,dim=['lat','lon'])
            
            bias = val - obs
            
            # Calculate RMSE error weighted by gridarea
            rmse = np.sqrt(masked_average(bias**2,weights=weights,dim=['lat','lon']))
    
            if contour:
                im = ax.contourf(bias['lon'],bias['lat'],bias,**kwargs)
                if lvl_bool:
                    print('SPECIFY LEVELS OR SUBPLOTS MAY NOT SHARE CONTOURS')
                    print('Levels: ', im.levels)
            else:
                im = bias.plot(ax=ax, **kwargs)
#             ax.set_title('Global Error: %.0f. RMS Error: %.0f' % (val_avg-obs_avg,rmse),fontsize=10)
            
            net_label = '%s Global Error: %.0f. RMS Error: %.0f' % (bold_label,val_avg-obs_avg,rmse)
            ax.set_title(net_label,fontsize=12) # JKS
                    
        else:
            weights = add_weights(val)['cell_weight']
            val_avg = masked_average(val,weights=weights,dim=['lat','lon'])
            if contour:
                im = ax.contourf(val['lon'],val['lat'],val,**kwargs)
                if lvl_bool:
                    print('SPECIFY LEVELS OR SUBPLOTS MAY NOT SHARE CONTOURS')
                    print('Levels: ', im.levels)
                    
            else:
                im = val.plot(ax=ax,**kwargs) # robust good?

            net_label = '%s Global Average: %.0f' % (bold_label,val_avg)
            ax.set_title(net_label,fontsize=12) # JKS

        add_map_features(ax) # testing turning this off
        return ax, im
    
    def band_biases(self, var):
        '''
        Calculate the model bias with respect to observations over 10 degree(ish) latitude bands.
        '''
        self.__check_for_vars(var)
        obs_source, obs_label = self.__get_data_source(var)
        out_dict = {}
        
        # Preprocessing observational data:
        obs = obs_source[var].mean('time', skipna=True)
        obs_avg = []
        for band in self.lat_bounds:
            # calculate weighted average for goccp. Open bottom is arbitrary. 
            # Remember weird mask convention! dah!
            goccp_weights = obs_source['cell_weight']
            goccp_mask = np.bitwise_or(obs_source['lat']<=band[0], 
                                       obs_source['lat']>band[1])
            goccp_avg = masked_average(obs, dim=['lat','lon'], 
                                       weights=goccp_weights, mask=goccp_mask)
            obs_avg.append(goccp_avg.values)
        print("obs_avg: ", obs_avg)
        
        print("Latitude Band Error (Model - GOCCP) for %s." % var)
        # Iterate over the cases:
        for k in self.cases:
            _biases = []
            _run = self.cases[k]
            print(k)
            _da = _run.case_da
            _val = _da[var].mean('time', skipna=True)
            _weights = _da['cell_weight']
            for band,_goccp in zip(self.lat_bounds,obs_avg):
                _mask = np.bitwise_or(_da['lat']<=band[0], 
                                      _da['lat']>band[1])
                _avg = masked_average(_val, dim=['lat','lon'], 
                                      weights=_weights, mask=_mask)
#                 print("Model: ", _avg.values, "Obs.: ", _goccp.value)
                print("%s error: %s" % (_run.label, (_avg - _goccp).values))
                _biases.append(_avg.values)
            out_dict[k] = _biases
#             print(_avg.values, _da['lat'])
        return out_dict

    def band_bias(self, var, band_bounds, season=None):
        '''
        Calculate the model bias with respect to observations over a specified latitude band.
        '''
        self.__check_for_vars(var)
        obs_source, obs_label = self.__get_data_source(var)
        out_dict = {}
        out_list = []
        # Preprocessing observational data:
        if season:
            obs = season_mean(obs_source[var])[self.seas_dict[season]]#.to_dataset(name=var)
        else:
            obs = obs_source[var].mean('time', skipna=True)
            # calculate weighted average for goccp. Open bottom is arbitrary. 
            # Remember weird mask convention! dah!
        goccp_weights = obs_source['cell_weight']
        goccp_mask = np.bitwise_or(obs_source['lat']<=band_bounds[0], 
                                       obs_source['lat']>band_bounds[1])
        goccp_avg = masked_average(obs, dim=['lat','lon'], 
                                    weights=goccp_weights, mask=goccp_mask)
        
#         print("Latitude Band Error (Model - GOCCP) for %s over %.2f - %.2f lat." % 
#               (var, band_bounds[0],band_bounds[1]))
        # Iterate over the cases:
        for k in self.cases:
            _biases = []
            _run = self.cases[k]
#             print(k)
            if season:
                _season_da = season_mean(_run.case_da[var])
                _da = _season_da[self.seas_dict[season]].to_dataset(name=var)
                _da = add_weights(_da)
                _val = _da[var]
            else:
                _da = _run.case_da
                _val = _da[var].mean('time', skipna=True)
            _weights = _da['cell_weight']
            _mask = np.bitwise_or(_da['lat']<=band_bounds[0], 
                                  _da['lat']>band_bounds[1])
            _avg = masked_average(_val, dim=['lat','lon'], 
                                  weights=_weights, mask=_mask)
#             print("%s error: %s" % (_run.label, (_avg - goccp_avg).values))
            out_list.append([_run.label, (_avg - goccp_avg).values])
        return out_list
        
    def plot_slf_isos(self):
        isos_plot = plt.figure(figsize=[10,10])
        plt.gca().invert_yaxis()

        # Plot satellite phase retrievals
        ic_caliop_weight = self.incloud_caliop_slf['cell_weight']
        ic_caliop_mask = np.bitwise_or(self.incloud_caliop_slf['lat']<66.667, 
                                       self.incloud_caliop_slf['lat']>82)
        ic_caliop_slf = 100*masked_average(self.incloud_caliop_slf['SLF'], dim=['lat','lon'],
                                           weights=ic_caliop_weight, mask=ic_caliop_mask)
        ic_caliop_stdev = 100*np.std(self.incloud_caliop_slf['SLF'].sel(lat=slice(66.667,90)), axis=(0,1))
        plt.errorbar(ic_caliop_slf, ic_caliop_slf['isotherm'], xerr=ic_caliop_stdev, 
                     label='CALIOP SLF IC', fmt='o-', color = 'black', zorder=0)
        
        
        caliop_weight = self.ct_caliop_slf['cell_weight']
        caliop_mask = np.bitwise_or(self.ct_caliop_slf['lat']<66.667, self.ct_caliop_slf['lat']>82)
        caliop_slf = 100*masked_average(self.ct_caliop_slf['SLF'], dim=['lat','lon'], 
                                        weights=caliop_weight, mask=caliop_mask)
        caliop_stdev = 100*np.std(self.ct_caliop_slf['SLF'].sel(lat=slice(66.667,90)), axis=(0,1))
        _line = plt.errorbar(caliop_slf, caliop_slf['isotherm'], xerr=caliop_stdev, label='CALIOP SLF',
                             fmt='o', color = 'black', zorder=0, linestyle='-', marker='D')

        labels = ['CALIOP']
        lines = [_line]
        for i,color in zip(self.cases, self.colors):
            _run = self.cases[i]
            _case = _run.case_da
            _weight = _case['cell_weight']#*_case['CT_CLD_ISOTM'] #Not sure about this weighting
            _mask = np.bitwise_or(_case['lat']<66.667, _case['lat']>82)
            
            # Cloudtop SLF part
            _case['CT_SLF_TAVG'] = _case['CT_SLF'].mean(dim = 'time', skipna=True)
            slf_ct = 100*masked_average(_case['CT_SLF_TAVG'], dim=['lat','lon'], weights=_weight, mask=_mask)
            err = np.array(slf_ct) - np.array(caliop_slf)
            rms_ct = np.sqrt(np.mean(np.square(err)))

            _line = plt.scatter(slf_ct, slf_ct['isotherms_mpc'] - 273.15, 
                                label=(_run.label+' RMSE: %.2f' % rms_ct), color=color, marker='D')
            
            # Bulk SLF part
            _case['SLF_ISOTM'] = (_case['SLFXCLD_ISOTM'] / _case['CLD_ISOTM'])
            slf_bulk = 100*masked_average(_case['SLF_ISOTM'], dim=['lat','lon', 'time'], weights=_weight, mask=_mask)
            err = np.array(slf_bulk) - np.array(ic_caliop_slf)       
            rms_bulk = np.sqrt(np.mean(np.square(err)))
            
            plt.scatter(slf_bulk, slf_bulk['isotherms_mpc'] - 273.15, 
                        label=(_run.label+' RMSE: %.2f' % rms_bulk), color=color)
#             labels.append(_run.label+' CT_RMSE: %.2f, Bulk_RMSE: %.2f' % (rms_ct, rms_bulk))
            labels.append(_run.label) #+' CT_RMSE: %.2f, Bulk_RMSE: %.2f' % (rms_ct, rms_bulk))
            lines.append(_line)
            
        plt.xlabel('SLF (%)', fontsize=18)
        plt.ylabel('Isotherm (C)', fontsize=18)
        plt.legend(lines, labels)
        plt.title('SLF trends with microphysical modifications', fontsize=24)


    def cloud_polar_plot(self, lat_range=[66,82]):
        '''
        This will create a polar plot to show the monthly evolution 
        of the cloud phase breakdown. In development.
        '''
        self.__check_for_vars(['CLDTOT_CAL','CLDTOT_CAL_LIQ','CLDTOT_CAL_ICE','CLDTOT_CAL_UN'])
        fig = plt.figure(figsize=[6*(len(self.cases)),10],dpi=200)
        # Rotate and reverse so that they read clockwise from the top
        months_dq = deque(self.months) 
        months_dq.rotate(-4) 
        months_dq = list(months_dq)
        months_dq.reverse()
        
        # Process CALIOP_GOCCP
        obs_da = self.goccp_data.sel(lat=slice(lat_range[0],lat_range[1]))
        _tot_vals = self.__average_and_wrap(obs_da['CLDTOT_CAL'])
        _liq_vals = self.__average_and_wrap(obs_da['CLDTOT_CAL_LIQ'])
        _ice_vals = self.__average_and_wrap(obs_da['CLDTOT_CAL_ICE'])
        _un_vals = self.__average_and_wrap(obs_da['CLDTOT_CAL_UN'])
        
        # top is pi/2, proceed counterclockwise
        theta = np.linspace(np.pi/2, -3/2 * np.pi, len(self.months)+1) # rotated and reversed
        for i,k in enumerate(self.cases):
            ax = fig.add_subplot(1,len(self.cases),i+1, projection='polar')
        
            lines, labels = plt.thetagrids(range(0, 360, int(360/len(months_dq))), (months_dq)) 

            # Plot observations from CALIOP-GOCCP
            obs_tot = plt.plot(theta, _tot_vals, marker='D', linestyle='dashed', color='orange')
            obs_liq = plt.plot(theta, _liq_vals, marker='D',linestyle='dashed',color='blue')
            obs_lice = plt.plot(theta, _liq_vals+_ice_vals, marker='D',linestyle='dashed',color='green')
            blank=np.empty(theta.shape)
            blank[:] = np.NaN
            plt.plot(theta,blank,marker='D',color='gray',label='CALIOP-GOCCP')

            # Plot model data from the current run
            _run = self.cases[k]
            _da = _run.case_da.sel(lat=slice(lat_range[0],lat_range[1]))
            tot_vals = self.__average_and_wrap(_da['CLDTOT_CAL'])
            liq_vals = self.__average_and_wrap(_da['CLDTOT_CAL_LIQ'])
            ice_vals = self.__average_and_wrap(_da['CLDTOT_CAL_ICE'])
            un_vals = self.__average_and_wrap(_da['CLDTOT_CAL_UN'])
#             tot_vals = np.roll(self.__average_and_wrap(_da['CLDTOT_CAL']),2)

            plt.plot(theta, liq_vals,color='blue',linewidth=0.5)
            plt.fill(theta, liq_vals, color='blue', alpha = 0.1, label="Liquid Cloud") #,hatch='//') 
            plt.plot(theta, liq_vals+ice_vals,color='green',linewidth=0.5) 
            plt.fill_between(theta, liq_vals, liq_vals+ice_vals,label="Ice Cloud",alpha=0.1,color='green')#, 'b', alpha = 0.1, label="Ice Cloud") 
            plt.plot(theta, tot_vals,color='orange',linewidth=0.5)
            plt.fill_between(theta, tot_vals, liq_vals+ice_vals,label="Undefined Phase",alpha=0.1,color='orange')#, 'b', alpha = 0.1, label="Ice Cloud",linewidth=0.5) 

            plt.title(_run.label, fontsize=15)    
            ax.set_rmax(100)
        plt.legend(loc='lower right', bbox_to_anchor=(0.8, -0.1, 0.5, 0.5))
        return fig
    
    def plot_cloud_sum(self, bias=False, season=None):
        '''
        Latitude plot of cloud phase components.
        '''
        self.__check_for_vars(['CLDTOT_CAL','CLDTOT_CAL_LIQ','CLDTOT_CAL_ICE'])
        fig, axes = plt.subplots(nrows=1,ncols=len(self.cases),figsize=[15,4])

        obs_da = self.goccp_data
        obs_handle = mlines.Line2D([], [], linestyle='dashed', color='gray', label='CALIOP-GOCCP')
        if season:
            obs_cldtot = self.__seasonalize_for_sum(obs_da, 'CLDTOT_CAL', season)
            obs_cldliq = self.__seasonalize_for_sum(obs_da, 'CLDTOT_CAL_LIQ', season)
            obs_cldice = self.__seasonalize_for_sum(obs_da, 'CLDTOT_CAL_ICE', season)
            
            handles = []
            
        else:
            obs_cldtot = obs_da['CLDTOT_CAL'].mean(['time','lon'])
            obs_cldliq = obs_da['CLDTOT_CAL_LIQ'].mean(['time','lon'])
            obs_cldice = obs_da['CLDTOT_CAL_ICE'].mean(['time','lon'])

        for k,ax in zip(self.cases,axes):
            _run = self.cases[k]
            _da = _run.case_da
            
            if season:
                _cldtot = self.__seasonalize_for_sum(_da, 'CLDTOT_CAL', season)
                _cldliq = self.__seasonalize_for_sum(_da, 'CLDTOT_CAL_LIQ', season)
                _cldice = self.__seasonalize_for_sum(_da, 'CLDTOT_CAL_ICE', season)
            else:
                _cldtot = _da['CLDTOT_CAL'].mean(['time','lon'])
                _cldliq = _da['CLDTOT_CAL_LIQ'].mean(['time','lon'])
                _cldice = _da['CLDTOT_CAL_ICE'].mean(['time','lon'])

            # Plot model output
            ax.plot(_cldliq['lat'], _cldliq, color='blue',linewidth=0.5)
            liq_fill = ax.fill_between(_cldliq['lat'], _cldliq, np.zeros(_cldliq.shape), color='blue', alpha = 0.1, label="Liquid Cloud")
            ax.plot(_cldice['lat'], _cldice+_cldliq,color='green',linewidth=0.5) 
            lice_fill = ax.fill_between(_cldliq['lat'], _cldliq, _cldice+_cldliq, label="Ice Cloud",alpha=0.1,color='green')
            ax.plot(_cldtot['lat'], _cldtot,color='orange',linewidth=0.5) 
            un_fill = ax.fill_between(_cldliq['lat'], _cldice+_cldliq, _cldtot, label="Undefined Phase",alpha=0.1,color='orange')

            # Plot observations
            obs_tot = ax.plot(obs_cldtot['lat'], obs_cldtot, linestyle='dashed', color='orange')
            obs_liq = ax.plot(obs_cldliq['lat'], obs_cldliq, linestyle='dashed',color='blue')
            obs_lice = ax.plot(obs_cldliq['lat'], obs_cldliq+obs_cldice, linestyle='dashed',color='green')
            blank=np.empty(obs_cldliq['lat'].shape)
            blank[:] = np.NaN
            obs_fill = ax.plot(obs_cldtot['lat'], blank, linestyle='dashed', color='gray', label='CALIOP-GOCCP')

        handles = [obs_handle, liq_fill, lice_fill, un_fill]

        ylabels = ['Cloud Cover (%)']; yax = [axes[0]]
        xlabels = self.case_labels; xax = axes
        self.add_labels(xlabels=xlabels, xaxes=xax, ylabels=ylabels, yaxes=yax)
        
#         plt.legend(handles=handles)
        plt.legend(handles=handles,loc='lower right', bbox_to_anchor=(0.8, -0.1, 0.5, 0.5))

    def plot_months_line(self, var, lat_range=[66,82], bias=False, ax=None, **kwargs):
            
        self.__check_for_vars(var)
        if ax:
            axes = ax
            plt.sca(ax)
            fig=None
        else:
            fig, axes = plt.subplots(nrows=1,ncols=1,figsize=[15,10])
        
        labels = []         
        obs_source, obs_label = self.__get_data_source(var)
        obs_source = obs_source.sel(lat=slice(lat_range[0],lat_range[1]))
        obs_vals = self.__average_and_wrap(obs_source[var],wrap=False)
        if not bias:
            obs_vals.plot(ax=axes,label=obs_label, color='black')
            labels.append(obs_label)
        
        lines = []
#         labels = []
        for k,color in zip(self.cases,self.colors):
            _run = self.cases[k]
            _da = _run.case_da.sel(lat=slice(lat_range[0],lat_range[1]))
            mon_vals = self.__average_and_wrap(_da[var],wrap=False)
            if bias:
                mon_vals = mon_vals - obs_vals
            _ln = mon_vals.plot(ax=axes,label=_run.label, color=color, **kwargs)
            _lbl = _run.label
            
#             lines.append(_ln)
            labels.append(_lbl)
            
        plt.xticks(np.arange(1,len(self.months)+1,1), self.months)
        plt.legend(labels) # JKS
        
        return fig
        
    def add_labels(self, xlabels=None, ylabels=None, xaxes=None, yaxes=None, height=1.05):
        '''
        Labelling function for 2d subplot grids.
        '''
        try:
            for ax,label in zip(xaxes,xlabels):
                ax.text(0.5, height, label.replace('_',' '), va='bottom', ha='center',
                        rotation='horizontal', rotation_mode='anchor',
                        transform=ax.transAxes, fontsize=14)
        except: pass

        try:
            for ax,label in zip(yaxes,ylabels):
                ax.text(-0.1, 0.5, label.replace('_',' '), va='bottom', ha='center',
                    rotation='vertical', rotation_mode='anchor',
                    transform=ax.transAxes, fontsize=14)
        except: pass
        
    def share_ylims(self, _axes):
        '''
        For 1D plots. Finds the global max and min so plots share bounds and are easier 
        to interpret.
        '''
        try:
            axes = _axes.flat # so single iteration works
        except:
            axes = _axes
        
        ymin, ymax = axes[0].get_ylim() # initialize values
        for ax in axes[1:]:
            _ymin, _ymax = ax.get_ylim()
            if _ymin < ymin: 
                ymin = _ymin
            if _ymax > ymax: 
                ymax = _ymax
                
        for ax in axes:
            ax.set_ylim([ymin,ymax])
            
        
    def share_clims(self, figure):
        '''
        For 2D plots. Finds the global max and min so plots share a colorbar 
        and are easier to interpret. In development.
        '''
        geoaxes = figure.axes

        qms = [] # to store QuadMesh object
        for i in geoaxes: # iterate over axes and find QuadMesh objects
    #         if isinstance(i, cpy.mpl.geoaxes.GeoAxesSubplot):
    # #             print(i.get_children())
            for j in i.get_children():
                # handle contour plots now
#                 print(type(j))
                if isinstance(j, mpl.collections.QuadMesh): #or isinstance(j, mpl.collections.PathCollection): 
                    qms.append(j)

        # Calculate global min and max values
        min,max = qms[0].get_clim() # initialize
        for _qm in qms:
            _clim = _qm.get_clim()
            if _clim[0] < min:
                min = _clim[0]
            if _clim[1] > max:
                max = _clim[1]
        
        # Extend by 10% (testing) JKS
        min = 1.1*min
        max = 1.1*max
                
        for _qm in qms:
            _qm.set_clim((min, max))
    
    def __get_data_source(self, var):
        ''' Figure out if a GOCCP or CERES-EBAF variable. Return correct dataset.'''
        if var in self.goccp_data.data_vars: # what is this
            data_source = self.goccp_data
            label = "CALIPSO" # was "GOCCP"
        elif var in self.ceres_data.data_vars:
            data_source = self.ceres_data
            label = "CERES-EBAF"
        elif var in self.isccp_data.data_vars:
            data_source = self.isccp_data
            label = "ISCCP"     
        elif var in self.misr_data.data_vars:
            data_source = self.misr_data
            label = "MISR"
        elif var in self.modis_data.data_vars:
            data_source = self.modis_data
            label = "MODIS"
        else:
            print("Could not find variable in GOCCP or CERES-EBAF datasets.")
            sys.exit() # bad?
        
        return data_source, label
            
    def __average_and_wrap(self,da,wrap=True):
        '''
        Helper function for cloud_polar_plot
        '''
        _dat = da.groupby('time.month').mean() # Create monthly averages
#         _dat = _dat.mean(['lat','lon'])#.values # Get np.array average, ! JKS use masked average
        
        # Gridcell weighting fix
        dat2 = add_weights(_dat)
        _weights = dat2['cell_weight']
            
        out_dat = masked_average(dat2, dim=['lat','lon'], weights=_weights)
        
        if wrap:
#             _dat = np.append(_dat, _dat[0]) # wrap
            out_dat = np.append(out_dat, out_dat[0]) # wrap
        
        return out_dat
#         return _dat
    
    def __seasonalize_for_sum(self, da, var, season):
        ''' Helper function for plot_cloud_sum'''
        season_da = season_mean(da[var])#.where(np.absolute(da[var]['lat'])<82)
        var_da = season_da[self.seas_dict[season]].to_dataset(name=var) # select correct season and revert to da
        out = var_da[var].mean('lon')
            
        return out
    
    def __check_for_vars(self,varlist):
        '''Check if each case has the requested variables. Add them if they are not present.''' 
        if not isinstance(varlist, list): # convert varlist to iterable if it isn't already
            varlist = [varlist]
                        
        for k in self.cases:
            _case = self.cases[k]
            for var in varlist:
                if not _case.case_da:
#                     print('Adding first variable %s' % var)
                    _case.add_tseries_var(var)
                if not var in _case.case_da.data_vars: # variable has not been loaded
#                     print('Adding variable %s' % var)
                    _case.add_tseries_var(var)
    
    def load_vars(self,varlist):
        '''Wrapper for internal function __check_for_vars'''
        self.__check_for_vars(varlist)

        
class Model_case:
    '''
    Class for models within a cloud metric class.
    
    '''
    
    def __init__(self, casedir, case, t_range=None):
        self.case_dir = casedir # General directory where all cases are stored
        self.case = case # The case name
        self.t_range = t_range
        
        self.add_ds()
        
        # Parse name to retrieve encoded values (dislike)
        try:
            _temp_parsed = self.case.split('_') # Parse string to get model parameters
            self.date, self.time, self.paramfile = _temp_parsed[:3]
            self.wbf_mult = np.float(_temp_parsed[-3]); self.inp_mult = np.float(_temp_parsed[-1])
            self.label = 'WBF: %.3f INP: %.3f' % (self.wbf_mult, self.inp_mult)
        
        except: # Handle non-standard case name formats
            self.label = self.case
            self.date, self.time, self.paramfile, self.wbf_mult, self.inp_mult = ['None','None','None','None','None']
        
    def add_ds(self): # try removing tseries here
        '''
        Load the CAM output.
        First try to identify a processed timeseries directory so that data can be loaded incrementally.
        If no processed data is found, use open_mfdataset (slow,bad).
        
        '''
        
        if os.path.exists("%s/%s/atm/proc/tseries/month_1/" % (self.case_dir,self.case)): # check for "processed tseries folder"
            print("Processed timeseries directory found for %s. Will load data as required." % self.case)
            self.tseries_input = True
            self.tseries_path = "%s/%s/atm/proc/tseries/month_1/" % (self.case_dir,self.case)
            self.alt_path = "%s/%s" % ('/glade/u/home/jonahshaw/w/archive/taylor_files/',self.case)# dirty hard-coded in for now 06/10/2021
            self.case_da = None
        else:
            try:
                print('Trying to load concatenated file for %s' % self.case)
                _ds = xr.open_dataset('%s%s/%s.nc' % (self.case_dir, self.case, self.case))
            except: # this probably wouldn't work anyway, at least on NIRD
                print('Failed, using xr.open_mfdataset.')
                self.geth0s() # create list of appropriate output files self.files
                _ds = xr.open_mfdataset(self.files)#, combine='by_coords') #, chunks={'lat':10}) #chunk?
            # Shift time-bounds to fix comparison with satellites (from Jen Kay)
            try:
                _ds['time'] = _ds['time_bnds'].isel(bnds=0)
            except:
                _ds['time'] = _ds['time_bnds'].isel(nbnd=0)

            if self.t_range:
                try:
                    ds = _ds.sel(time=slice(*self.t_range)) # flexible date selection with fun unpack syntax
                except:
                    print("Year format %s - %s not found." % (*self.t_range,))
            ds = add_weights(ds)

            self.case_da = ds
            self.tseries_input = False
        
        print("%s load successfully." % self.case)
        
    def geth0s(self):
        prepath = '%s%s/atm/hist/' % (self.case_dir, self.case)
        allfiles = os.listdir(prepath)
        self.files = [prepath + i for i in allfiles if "h0" in i]
        
    def add_tseries_var(self, var):
        '''
        Function adds timeseries data to the case_da attribute, creating
        '''
        # Get file path for variable timeseries
        tseries_files = os.listdir(self.tseries_path)
#         key_str = ".cam.h0.%s." % var # long enough to hopefully make this a unique identifier
        key_str = ".%s." % var # long enough to hopefully make this a unique identifier
        varfiles = [self.tseries_path + x for x in tseries_files if key_str in x]
        if len(varfiles) != 1:
            print('Not able to find unique timeseries file for %s in %s' % (var, self.label))
            print('Looking in alternate path: %s' % (self.alt_path))
            alt_files = os.listdir(self.alt_path)
            varfiles = [self.alt_path + '/' + x for x in alt_files if key_str in x]
            if len(varfiles) != 1:
                print('Alternate path also failed.')
                return
            else:
                print('Found in the alternate path.')

        varfile = varfiles[0]
        var_da = xr.open_dataset(varfile)
        try:
            var_da['time'] = var_da['time_bnds'].isel(bnds=0)
        except:
            var_da['time'] = var_da['time_bnds'].isel(nbnd=0) # wut. try-try again?

        if self.t_range:
            try:
                var_da = var_da.sel(time=slice(*self.t_range))
            except:
                print("Year format %s - %s not found." % (*self.t_range,))
                
        if self.case_da is None:
            # create self.case_da with var
            self.case_da = var_da
        else:
#             self.case_da = xr.merge([self.case_da, xr.open_dataset(varfile)]) # JKS testing
            self.case_da = xr.merge([self.case_da, var_da],compat='override') # JKS testing
mpl.rc('text', usetex=False)