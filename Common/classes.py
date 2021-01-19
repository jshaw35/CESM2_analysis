from imports import *
from functions import *

np.seterr(divide='ignore', invalid='ignore') # this maybe isn't working??


class CT_SLF_Metric(object):
    '''
    Used by JShaw for visualizing SLF data along isotherms from cloudtop and in-cloud metrics
    Initialized with a primary directory for locating data files and optional number of timesteps to use.
    Could be updated for a more arbitrary selection of time period.
    
    '''
    
    def __init__(self, casedir, timesteps=None, sel_month=None):
        np.seterr(divide='ignore', invalid='ignore') # testing to remove error message
        for cls in reversed(self.__class__.mro()):
            if hasattr(cls, 'init'):
                cls.init(self, casedir, timesteps, sel_month)

    def init(self, casedir, time_steps, sel_month):
        # Catch directory issues from the get-go
        if not os.path.exists(casedir):
            print('The case directory %s does not exist' % casedir)
#            sys.exit() # Is this a problem, does it need to be here?
            
        self.case_dir = casedir
        self.origin = None
        self.cases = {}
        self.time_steps = time_steps
        self.month = sel_month
#         self.colors = ['red','orange','yellow','green','blue','purple','pink']
#         self.colors = ['purple','red','orange','yellow','green','blue','pink']
#         self.colors = ['blue','darkorange','green','orange','yellow','green','blue','pink']
        self.colors = sns.color_palette("colorblind")

        self.seas_dict = {"DJF":0,"JJA":1,"MAM":2,"SON":3}
        
        try:
            self.ct_caliop_slf = xr.open_dataset('caliop_olimpia_new/netcdf_format/ct_slfs_annual.nc')
#             self.ct_caliop_slf = xr.open_dataset('caliop_olimpia/ct_slf_olimpia/cloudtop_slfs.nc')
        except:
            print('Could not load cloudtop cloud CALIOP slfs from caliop_olimpia/ct_slf_olimpia/cloudtop_slfs.nc')
        try:
            self.incloud_caliop_slf = xr.open_dataset('caliop_olimpia_new/netcdf_format/bulk_slfs_annual.nc')
#             self.incloud_caliop_slf = xr.open_dataset('caliop_olimpia/incloud_slf_olimpia/incloud_slfs.nc')
        except: 
            print('Could not load incloud cloud CALIOP slfs from caliop_olimpia/incloud_slf_olimpia/incloud_slfs.nc')            
        try:
            self.seasonal_ct_slf = xr.open_dataset('caliop_olimpia_new/netcdf_format/ct_slfs_seasonal.nc')
#             self.seasonal_ct_slf = xr.open_dataset('caliop_olimpia/seasonal_data/cloudtop_slfs_seasonal.nc')
        except:
            print('Could not load seasonal cloudtop SLFs')
        try:
            self.seasonal_incloud_slf = xr.open_dataset('caliop_olimpia_new/netcdf_format/bulk_slfs_seasonal.nc')
#             self.seasonal_incloud_slf = xr.open_dataset('caliop_olimpia/seasonal_data/incloud_slfs_seasonal.nc')
        except:
            print('Could not load seasonal incloud SLFs')
            
    def set_origin(self, case): # case is a string of the case name
        if case not in self.cases: # check if it is already in the dictionary
            self.add_case(case) 
        self.origin = self.cases[case] # should have a check for the object type
        # calculate values and set them to variables here
        
    def add_case(self, case, path=None, label=None): # Dictionary can be overwritten, should be fine
        if path == None:
            self.cases[case] = CT_SLF_case(self.case_dir, case, timesteps=self.time_steps, months=self.month, label=label) # Add a ct_slf_case object to the cases dictionary
        else:
            self.cases[case] = CT_SLF_case(path, case, timesteps=self.time_steps, months=self.month, label=label) # Add a ct_slf_case object to the cases dictionary
    
    def get_case(self, case):
        if case not in self.cases: # check if it is already in the dictionary
            print('Could not find %s in the SLF_Metric object' % case)
            return None
        else:
            return self.cases[case] # should have a check for the object type
    
    def get_cases(self):
        return self.cases
#        return list(self.cases.keys()) # This is a dict_keys object and should be an iterable
    
    def get_origin(self):
        return self.origin
    
    def plot_vectors(self):
        vector_plot = plt.figure(figsize=[10,10])

        origin = [self.origin.error_slope, self.origin.error_20c] # origin point (x,y)    
        for i in self.cases:
            _case = self.cases[i]
            # xy is the endpoint, xytext is the origin
            label = 'WBF: %.3f INP: %.3f, %.3f SLF RMS error' % (_case.wbf_mult, _case.inp_mult, _case.rms_error)
            plt.scatter(_case.error_slope, _case.error_20c, label=label)
            plt.annotate('', xy=(_case.error_slope, _case.error_20c), xytext=(origin[0], origin[1]), 
                        arrowprops=dict(facecolor='black', shrink=0., width=0.5, lw=1.0),
                        label = label)

        plt.scatter([[0]], [[0]], color = 'black')
        plt.xlabel('Error Slope (SLF % / C)')
        plt.ylabel('0C Isotherm SLF Error')
        plt.grid()
        plt.legend()

        return vector_plot
    
    def plot_isos_bulk(self):
        isos_plot = plt.figure(figsize=[10,10])
        plt.gca().invert_yaxis()

        # use old observational data
       # plt.plot(self.bulk_caliop_slf['CALIOP_70S-90S'], self.bulk_caliop_slf['Isotherm'], 
        #         color='black', label='CALIOP', linestyle='-', marker='D')

        caliop_weight = self.incloud_caliop_slf['cell_weight']
        caliop_mask = np.bitwise_or(self.incloud_caliop_slf['lat']<70, self.incloud_caliop_slf['lat']>82)
        caliop_slf = 100*masked_average(self.incloud_caliop_slf['SLF'], dim=['lat','lon'], weights=caliop_weight, mask=caliop_mask)
        caliop_stdev = 100*np.std(self.incloud_caliop_slf['SLF'].sel(lat=slice(70,90)), axis=(0,1))
        plt.errorbar(caliop_slf, caliop_slf['isotherm'], xerr=caliop_stdev, label='CALIOP SLF', fmt='o-', marker='D', color = 'black', zorder=0)
        
        for i in self.cases:
            _run = self.cases[i]
            _case = _run.case_da
            _case['SLF_ISOTM'] = (_case['SLFXCLD_ISOTM'] / _case['CLD_ISOTM'])
            try:
                _case['SLF_ISOTM'] = (_case['SLFXCLD_ISOTM'] / _case['CLD_ISOTM']).sel(time=slice('0001-04-01',
                                    '0002-03-01'))
            except:
                _case['SLF_ISOTM'] = (_case['SLFXCLD_ISOTM'] / _case['CLD_ISOTM']).sel(time=slice('2000-04-01',
                                    '2001-03-01'))
                
            weight = _case['cell_weight']
            mask = np.bitwise_or(_case['lat']<70, _case['lat']>82) # caliop doesn't go north of 82lat
            slf = 100*masked_average(_case['SLF_ISOTM'], dim=['lat','lon', 'time'], weights=weight, mask=mask)
            
            err = np.array(slf)[::2][:-1] - np.array(self.bulk_caliop_slf['CALIOP_70S-90S'][:-1])
            rms = np.sqrt(np.mean(np.square(err)))
            
            _case['SLF_ISOTM_TAVG'] = _case['SLF_ISOTM'].mean(dim = 'time', skipna=True)
            
            stdev = 100*np.std(_case['SLF_ISOTM_TAVG'].sel(lat=slice(70,90)), axis=(1,2))
            plt.scatter(slf, slf['isotherms_mpc'] - 273.15, label=(_run.label+' RMSE: %.2f' % rms))
            
        plt.xlabel('SLF Error (%)')
        plt.ylabel('Isotherm (C)')
        plt.legend()
        plt.title('Bulk cloud SLF trends with INP and WBF modifications', fontsize=16)
        
        return isos_plot
    
    def plot_isos_ct(self):
        isos_plot = plt.figure(figsize=[10,10])
        plt.gca().invert_yaxis()

        caliop_weight = self.ct_caliop_slf['cell_weight']
        caliop_mask = np.bitwise_or(self.ct_caliop_slf['lat']<66.667, self.ct_caliop_slf['lat']>82)
        caliop_slf = 100*masked_average(self.ct_caliop_slf['SLF'], dim=['lat','lon'], weights=caliop_weight, mask=caliop_mask)
        caliop_stdev = 100*np.std(self.ct_caliop_slf['SLF'].sel(lat=slice(66.667,90)), axis=(0,1))
        plt.errorbar(caliop_slf, caliop_slf['isotherm'], xerr=caliop_stdev, label='CALIOP SLF', fmt='o', marker='D', color = 'black', zorder=0)

        for j, i in enumerate(self.cases):
            _run = self.cases[i]
            _case = _run.case_da
            _case['CT_SLF_TAVG'] = _case['CT_SLF'].mean(dim = 'time', skipna=True)

            weight = _case['cell_weight']*_case['CT_CLD_ISOTM'] #Not sure about this weighting
            mask = np.bitwise_or(_case['lat']<70, _case['lat']>90)
            slf = 100*masked_average(_case['CT_SLF_TAVG'], dim=['lat','lon', 'time'], weights=weight, mask=mask)
            err = np.array(slf) - np.array(caliop_slf)
            rms = np.sqrt(np.mean(np.square(err)))

            plt.scatter(slf, slf['isotherms_mpc'] - 273.15, label=(_run.label+' RMSE: %.2f' % rms), zorder=j)
        
        plt.xlabel('SLF Error (%)')
        plt.ylabel('Isotherm (C)')
        plt.legend()
        plt.title('Cloudtop SLF trends with INP and WBF modifications', fontsize=16)
                               
        return isos_plot
    
    def plot_isos_all(self, lat_range=[66.667,82],season=None, **kwargs):        
#         isos_plot = plt.figure(figsize=[10,10])
#         print(kwargs['ax'])
        try:
#             kwargs['ax'].invert_yaxis()
            plt.sca(kwargs['ax'])
            isos_plot = None
        except: 
            isos_plot = plt.figure(figsize=[10,10])
        plt.gca().invert_yaxis()

        if season: # handle looking at specific seasons
            ic_season_da = season_mean(self.seasonal_incloud_slf['SLF'])
            ct_season_da = season_mean(self.seasonal_ct_slf['SLF'])
            _ic_slf = ic_season_da.sel(season=season).to_dataset(name='SLF')
            _ct_slf = ct_season_da.sel(season=season).to_dataset(name='SLF')
            ic_slf = add_weights(_ic_slf)
            ct_slf = add_weights(_ct_slf)
        else:
            ic_slf = self.incloud_caliop_slf
            ct_slf = self.ct_caliop_slf
        
        # Plot satellite phase retrievals
        ic_caliop_weight = ic_slf['cell_weight']
        ic_caliop_mask = np.bitwise_or(ic_slf['lat']<lat_range[0], 
                                       ic_slf['lat']>lat_range[1])
        ic_caliop_slf = 100*masked_average(ic_slf['SLF'], dim=['lat','lon'], weights=ic_caliop_weight, mask=ic_caliop_mask)
        ic_caliop_stdev = 100*np.std(ic_slf['SLF'].sel(
                                     lat=slice(lat_range[0],lat_range[1])),axis=(0,1))
        ic_line = plt.errorbar(ic_caliop_slf, ic_caliop_slf['isotherm'], xerr=ic_caliop_stdev, label='CALIOP SLF IC', fmt='o-', color='black', zorder=0, linestyle='dashed') # (0,(5,1)) '--'
        
        
        caliop_weight = ct_slf['cell_weight']
        caliop_mask = np.bitwise_or(ct_slf['lat']<lat_range[0], 
                                    ct_slf['lat']>lat_range[1])
        caliop_slf = 100*masked_average(ct_slf['SLF'], dim=['lat','lon'], weights=caliop_weight, mask=caliop_mask)
        caliop_stdev = 100*np.std(ct_slf['SLF'].sel(lat=slice(lat_range[0],lat_range[1])), axis=(0,1))
        _line = plt.errorbar(caliop_slf, caliop_slf['isotherm'], xerr=caliop_stdev, label='CALIOP SLF', color = 'black', zorder=0, linestyle='-', marker='D')# fmt='o',
               
        labels = ['CALIOP Cloudtop SLF','CALIOP Cloudbulk SLF']
        lines = [_line, ic_line] 
        for i,color in zip(self.cases, self.colors): # [3:] jks
            _run = self.cases[i]
            _case = _run.case_da
            
            # Cloudtop SLF part
            if season: # handle looking at specific seasons
                _case['CT_SLF_TAVG'] = season_mean(_case['CT_SLF']).sel(season=season)
            else:
                _case['CT_SLF_TAVG'] = _case['CT_SLF'].mean(dim = 'time', skipna=True)
            # this order of things must be thought about.
            weight_ct = _case['cell_weight']#*_case['CT_CLD_ISOTM'] #Not sure about this weighting
            mask = np.bitwise_or(_case['lat']<lat_range[0], _case['lat']>lat_range[1])
            slf_ct = 100*masked_average(_case['CT_SLF_TAVG'], dim=['lat','lon'], weights=weight_ct, mask=mask) # error here related to having 'time' or not
            err = np.array(slf_ct) - np.array(caliop_slf)
            rms_ct = np.sqrt(np.mean(np.square(err)))

            _line = plt.scatter(slf_ct, slf_ct['isotherms_mpc'] - 273.15, label=(_run.label+' RMSE: %.2f' % rms_ct), color=color, marker='D')
            _line2 = plt.plot(slf_ct, slf_ct['isotherms_mpc'] - 273.15, label=(_run.label+' RMSE: %.2f' % rms_ct), color=color, linestyle='-')
            
            # Bulk SLF part
            if season: # handle looking at specific seasons
                _case['SLF_ISOTM'] = season_mean(_case['SLFXCLD_ISOTM'] / _case['CLD_ISOTM']).sel(season=season)
            else:
                _case['SLF_ISOTM'] = (_case['SLFXCLD_ISOTM'] / _case['CLD_ISOTM']).mean(dim = 'time', skipna=True)
            weight_bulk = _case['cell_weight']
            slf_bulk = 100*masked_average(_case['SLF_ISOTM'], dim=['lat','lon'], weights=weight_bulk, mask=mask)
            err = np.array(slf_bulk) - np.array(ic_caliop_slf)       
            rms_bulk = np.sqrt(np.mean(np.square(err)))
            
            plt.scatter(slf_bulk, slf_bulk['isotherms_mpc'] - 273.15, label=(_run.label+' RMSE: %.2f' % rms_bulk), color=color)
            plt.plot(slf_bulk, slf_bulk['isotherms_mpc'] - 273.15, label=(_run.label+' RMSE: %.2f' % rms_bulk), color=color,linestyle='dashed')
            labels.append(_run.label+' (%.2f, %.2f)' % (rms_ct, rms_bulk))
#             labels.append(_run.label+' CT_RMSE: %.2f, Bulk_RMSE: %.2f' % (rms_ct, rms_bulk))
            lines.append(_line)
        plt.xlim([-5,105])    
        plt.xlabel('SLF (%)', fontsize=18)
        plt.ylabel('Isotherm (C)', fontsize=18)
        plt.legend(lines, labels, loc='lower left',fontsize=15) # jks
        plt.title('SLF trends with microphysical modifications', fontsize=24)
        
        return isos_plot
    
    def plot_single_var(self, var_str, lat_range=[66.667,90]):
        '''
        Plot a single variable against 'lev' for loaded cases. Argument is the variable string name.
        '''
        
        test_plot = plt.figure(figsize=[10,10])
        
        test_plot.gca().invert_yaxis()
        for _casename, _col in zip(self.cases, self.colors):
            _case = self.cases[_casename]
            _da = _case.case_da
#            print(_da.variables)
            weight = _da['cell_weight']
            working_da = _da[var_str].sel(lat=slice(lat_range[0],lat_range[1])) # select variable and lat range
            working_da = working_da.mean(dim='time') # select time or average by it
#            weight = working_da['cell_weight']
            avg_var = masked_average(working_da, dim=['lat','lon'],weights=weight) # finally actually compute on the bounded variable of interest
            plt.plot(avg_var, avg_var['lev'], label=_case.label + " " + _case.date, color=_col)
        
#        return _da
        xaxis_label = '%s (%s)' % (var_str, _da[var_str].units)

        plt.legend(loc='upper right')
        test_plot.suptitle(_da[var_str].long_name, fontsize=24)
        plt.xlabel(xaxis_label); plt.ylabel('Hybrid Sigma pressure (hPa)');

        return test_plot

    def plot_var_ratio(self, num_var, denom_var, lat_range=[66.667,90]):
        '''
        Plot the ratio of two variables. This could be somehow combined with plot_single_var
        '''
        
        test_plot = plt.figure(figsize=[10,10])
        
        test_plot.gca().invert_yaxis()
        for _casename, _col in zip(self.cases, self.colors):
            _case = self.cases[_casename]
            _da = _case.case_da
#            print(_da.variables)
            weight = _da['cell_weight']
            working_da = (_da[num_var] / _da[denom_var]).sel(lat=slice(lat_range[0],
                                                                       lat_range[1])) # select variable and lat range
            working_da = working_da.mean(dim='time') # select time or average by it
            avg_var = masked_average(working_da, dim=['lat','lon'],weights=weight) # finally actually compute on the bounded variable of interest
            plt.plot(avg_var, avg_var['lev'], label=_case.label + " " + _case.date, color=_col)
        
#        return _da
        xaxis_label = '%s / %s' % (num_var, denom_var) # units issue here

        plt.legend(loc='upper right')
        test_plot.suptitle(_da[num_var].long_name, fontsize=24)
        plt.xlabel(xaxis_label, fontsize=16); plt.ylabel('Hybrid Sigma pressure (hPa)', fontsize=16);

        return test_plot
    
    
    def plot_parameterspace(self):
        parameterspace_plot = plt.figure(figsize=[10,10])       
#        origin = [test.get_origin().wbf_mult, test.get_origin().inp_mult] # origin point (x,y)

        plt.yscale('log')
        plt.xscale('log')

        for i in self.cases:
            _case = self.cases[i]
            label = 'WBF: %.3f INP: %.3f, %.3f SLF RMS error' % (_case.wbf_mult, _case.inp_mult, _case.rms_error)
            
#            if _case.inp_mult == 0: Not sure how to catch this error
            plt.annotate('', xy=(_case.wbf_mult, _case.inp_mult), xytext=(self.origin.wbf_mult, self.origin.inp_mult), 
                        arrowprops=dict(facecolor='black', shrink=0., width=0.5, lw=1.0),
                        )
            plt.scatter(_case.wbf_mult, _case.inp_mult, label=label)

        plt.scatter([self.origin.wbf_mult], [self.origin.inp_mult], label = 'Initial Parametrization')
        plt.scatter([[0]], [[0]], color = 'black')
        #plt.xlim((1e-3,1e2));
        #plt.ylim((1e-3,1e2))
        plt.xlabel('WBF Multiplier')
        plt.ylabel('INP Multiplier')
        plt.grid()
#        plt.legend()
        
        return parameterspace_plot

    def dual_parameterspace(self):
        parameterspace_plot = plt.figure(figsize=[10,10])       
        
        # Satellite values
        ic_caliop_weight = self.incloud_caliop_slf['cell_weight']
        ic_caliop_mask = np.bitwise_or(self.incloud_caliop_slf['lat']<70, self.incloud_caliop_slf['lat']>82)
        ic_caliop_slf = 100*masked_average(self.incloud_caliop_slf['SLF'], dim=['lat','lon'], weights=ic_caliop_weight, mask=ic_caliop_mask)
        
        caliop_weight = self.ct_caliop_slf['cell_weight']
        caliop_mask = np.bitwise_or(self.ct_caliop_slf['lat']<66.667, self.ct_caliop_slf['lat']>82)
        caliop_slf = 100*masked_average(self.ct_caliop_slf['SLF'], dim=['lat','lon'], weights=caliop_weight, mask=caliop_mask)
        
        for i, color in zip(self.cases, self.colors):
            _run = self.cases[i]
            _case = _run.case_da
            
            _case['CT_SLF_TAVG'] = _case['CT_SLF'].mean(dim = 'time', skipna=True)
            # this order of things must be thought about.
            weight = _case['cell_weight']#*_case['CT_CLD_ISOTM'] #Not sure about this weighting
            mask = np.bitwise_or(_case['lat']<70, _case['lat']>82)
            slf_ct = 100*masked_average(_case['CT_SLF_TAVG'], dim=['lat','lon'], weights=weight, mask=mask)
            
            _case['SLF_ISOTM'] = (_case['SLFXCLD_ISOTM'] / _case['CLD_ISOTM'])
            slf_bulk = 100*masked_average(_case['SLF_ISOTM'], dim=['lat','lon', 'time'], weights=weight, mask=mask)

            # Calculate errors
            err_ic = np.mean(np.array(slf_bulk) - np.array(ic_caliop_slf))
            err_ct = np.mean(np.array(slf_ct) - np.array(caliop_slf))
            
            change_ic = np.mean(np.array(slf_ct) - np.array(caliop_slf))
            change_ct = np.mean(np.array(slf_ct) - np.array(caliop_slf))
            
            print(err_ic, err_ct)
            
            plt.plot(err_ic, err_ct, marker='o', color=color, label=_run.label)
        
        plt.scatter(0,0,marker=(5, 1), color='gold',label='Goal',s=200)
        plt.ylabel('Avg. Error in Cloudtop SLF')
        plt.xlabel('Avg. Error in In-Cloud SLF')
        plt.grid(True)
        plt.legend()

    def parameterspace_change(self):
        parameterspace_plot = plt.figure(figsize=[10,10])       
        
        # Satellite values
        ic_caliop_weight = self.incloud_caliop_slf['cell_weight']
        ic_caliop_mask = np.bitwise_or(self.incloud_caliop_slf['lat']<70, self.incloud_caliop_slf['lat']>82)
        ic_caliop_slf = 100*masked_average(self.incloud_caliop_slf['SLF'], dim=['lat','lon'], weights=ic_caliop_weight, mask=ic_caliop_mask)
        
        caliop_weight = self.ct_caliop_slf['cell_weight']
        caliop_mask = np.bitwise_or(self.ct_caliop_slf['lat']<66.667, self.ct_caliop_slf['lat']>82)
        caliop_slf = 100*masked_average(self.ct_caliop_slf['SLF'], dim=['lat','lon'], weights=caliop_weight, mask=caliop_mask)
        
        oo = self.origin
        o_bulkslf, o_ctslf = self.arctic_slf(oo)
        
        #err_ic = np.mean(np.array(_bulkslf) - np.array(ic_caliop_slf))
        #err_ct = np.mean(np.array(_ctslf) - np.array(caliop_slf))
        
        for i, color in zip(self.cases, self.colors):
            _run = self.cases[i]
            
            _bulkslf, _ctslf = self.arctic_slf(_run)

            # Calculate errors
            #err_ic = np.mean(np.array(_bulkslf) - np.array(ic_caliop_slf))
            #err_ct = np.mean(np.array(_ctslf) - np.array(caliop_slf))
            
            change_ic = np.mean(np.array(_bulkslf) - np.array(o_bulkslf))
            change_ct = np.mean(np.array(_ctslf) - np.array(o_ctslf))
            
            #print(err_ic, err_ct)
            
            plt.plot(change_ic, change_ct, marker='o', color=color, label=_run.label)
        
        plt.ylabel('Avg. Change in Cloudtop SLF')
        plt.xlabel('Avg. Change in In-Cloud SLF')
        plt.grid(True)
        plt.legend()
        
        
    def arctic_slf(self, case):
        '''
        Function for use within the class to calculate SLF values in the CALIOP Arctic. Returns SLF by isotherm for both metrics.
        '''
        _case = case.case_da
            
        _case['CT_SLF_TAVG'] = _case['CT_SLF'].mean(dim = 'time', skipna=True)
        # this order of things must be thought about.
        weight = _case['cell_weight']#*_case['CT_CLD_ISOTM'] #Not sure about this weighting
        mask = np.bitwise_or(_case['lat']<70, _case['lat']>82)
        slf_ct = 100*masked_average(_case['CT_SLF_TAVG'], dim=['lat','lon'], weights=weight, mask=mask)
            
        _case['SLF_ISOTM'] = (_case['SLFXCLD_ISOTM'] / _case['CLD_ISOTM'])
        slf_bulk = 100*masked_average(_case['SLF_ISOTM'], dim=['lat','lon', 'time'], weights=weight, mask=mask)
        
        return slf_bulk, slf_ct
        
class CT_SLF_case:
    
    def __init__(self, casedir, case, timesteps, months, label):
        self.case_dir = casedir # General directory where all cases are stored
        self.case = case # The origin case name (wbf = 1, slf = 1)
        self.time_steps = timesteps
        self.months = np.array(months) # How does this work?
        
        self.add_ds()
        
        _temp_parsed = self.case.split('_') # Parse string to get model parameters
        try:
            self.date, self.time, self.paramfile = _temp_parsed[:3]
            self.wbf_mult = np.float(_temp_parsed[-3]); self.inp_mult = np.float(_temp_parsed[-1])
            self.label = 'WBF: %.3f INP: %.3f' % (self.wbf_mult, self.inp_mult)
        except: # Handle non-standard case name formats
            self.label = self.case
            self.date, self.time, self.paramfile, self.wbf_mult, self.inp_mult = ['None','None','None','None','None']
        if label: self.label=label
        
    def add_ds(self):
    # create list of appropriate files
        try:
            print('Trying to load concatenated file for %s' % self.case)
            _ds = xr.open_dataset('%s%s/%s.nc' % (self.case_dir, self.case, self.case))
        except:
            print('Failed, creating mfdataset.')
            self.geth0s()
            _ds = xr.open_mfdataset(self.files)#, combine='by_coords') #, chunks={'lat':10}) #chunk?
            
          #  print('Could not load from output files.')
                
#         try:
#             _ds = xr.open_dataset('%s%s/%s.nc' % (self.case_dir, self.case, self.case))
#         except:
#             print('Loading abridged dataset "slfvars" for %s.' % self.case)
#             _ds = xr.open_dataset('%s%s/%s_slfvars.nc' % (self.case_dir, self.case, self.case))
        # Do I only want to use this processed file? Probably variable-wise.
#        _ds = xr.open_mfdataset('%s%s/%s_slfvars.nc' % (self.case_dir, self.case, self.case), combine='by_coords')

        if self.time_steps != None: # if user has specified the number of steps to be used just index
            print("Just taking the first %s timesteps" % self.time_steps)
            ds = _ds.isel(time=slice(0,self.time_steps))
        
        elif self.months != None:  # select the desired months, 2 would mean all Februarys
            ds = _ds.sel(time=self.month_check(_ds['time.month']))
        
        elif (len(_ds['time']) > 12): # pick wind-up so months are equally weighted 20200428
            fake_windup = len(_ds['time']) % 12
            ds = _ds.sel(time=slice(_ds['time'][fake_windup],_ds['time'][-1]))
#            try:
#                ds = _ds.sel(time=slice('0001-04-01', '0002-03-01'))
#            except:
#                ds = _ds.sel(time=slice('2000-04-01', '2001-03-01'))                
        else: ds = _ds
            
        ds = add_weights(ds) # still has time here

        ds['CT_SLF'] = ds['CT_SLFXCLD_ISOTM']/ds['CT_CLD_ISOTM']
        ct_slf_noresm = ds['CT_SLF']

        ds['CT_SLF_ISOTM_AVG'] = ds['CT_SLF'].mean(dim = 'time', skipna=True)
        self.case_da = ds
        
    def month_check(self, month_data):
        _bools = []
        for i in month_data:
            try:
                _bools.append(i in self.months)
                
            except: pass
        
        return _bools
    
    def geth0s(self):
        prepath = '%s%s/atm/hist/' % (self.case_dir, self.case)
        allfiles = os.listdir(prepath)
        self.files = [prepath + i for i in allfiles if "h0" in i]
        
        
class SatComp_Metric(object):
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
                return None
            try:
                self.__load_CERESEBAF()
            except:
                print('Failed to load CERES-EBAF data.')
                
    
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
                                         'CLDHGH_CAL_ICE']}
        self.seasons = ['DJF','MAM','JJA','SON']
        self.projdict = {'PlateCarree':ccrs.PlateCarree(),'Arctic':ccrs.NorthPolarStereo(),
                         'Mollweide':ccrs.Mollweide(), "Antarctic":ccrs.SouthPolarStereo()}
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
            'toa_sw_all_mon':'FSNT','toa_lw_all_mon':'FLNT', # FSNT is wrong (needs reverse?)
            'toa_sw_clr_t_mon':'FSNTC','toa_lw_clr_t_mon':'FLNTC', # _t_ means that it is derived from all observations and is more consistent with the model variable
            'toa_cre_sw_mon':'SWCF','toa_cre_lw_mon':'LWCF',
            'sfc_lw_down_all_mon':'FLDS','sfc_sw_down_all_mon':'FSDS',
            'sfc_net_sw_all_mon':'FSNS','sfc_net_lw_all_mon':'FLNS',
            'sfc_net_lw_clr_t_mon':'FLNSC','sfc_net_sw_clr_t_mon':'FSNSC' # The longwave needs to be reversed
#             'toa_net_all_mon':'FNNTOA','toa_net_clr_c_mon':'FNNTOAC'
        }
        self.lat_bounds = [[-82,-70],[-70,-60],[-60,-50],[-50,-40],[-40,-30],
              [-30,-20],[-20,-10],[-10,0],[0,10],[10,20],[20,30],
              [30,40],[40,50],[50,60],[60,70],[70,82]]
        
        self.seas_dict = {"DJF":0,"JJA":1,"MAM":2,"SON":3}
        
        self.var_label_dict = {'CLDTOT_CAL':'Total Cloud \n Fraction (%)', 
                               'CLDTOT_CAL_LIQ':'Liquid Cloud \n Fraction (%)',
                               'CLDTOT_CAL_ICE':'Ice Cloud \n Fraction (%)'}
        
        self.months = ['J','F','M','A','M','J','J','A','S','O','N','D'] # month initials
        
        
        
    def __load_GOCCP_data(self):
        '''
        Load GOCCP data for model comparison. Relabel variables to match model output.
        
        Original data source:
        ftp://ftp.climserv.ipsl.polytechnique.fr/cfmip/GOCCP_v3/2D_Maps/grid_2x2xL40/
        '''
        print('Loading GOCCP data...', end = '')
        
        goccp_dir = 'GOCCP_data/2Ddata/f19_tn14_interpolation/'
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
        _goccp_data = add_weights(_goccp_data)

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
        _ceres_data = xr.open_dataset('CERES_EBAF/CERES_EBAF_SRFFLX_CRE_2X2.nc')
#         _ceres_data = xr.open_dataset('CERES_EBAF/CERES_EBAF_FLX_CRE_2X2.nc')
        _ceres_data = _ceres_data.sel(time=slice('2009-06-01', '2013-05-31'))
        _ceres_data = add_weights(_ceres_data)
        # Flip FLNSC so it matches model convention (net-LW is down, net-SW is up)
        _ceres_data['sfc_net_lw_clr_t_mon'] = -1*_ceres_data['sfc_net_lw_clr_t_mon']
        
        # Rename variables so they will match standard CAM variables names and index correctly
        self.ceres_data = _ceres_data.rename(self.ceres_var_dict)
        print('done.')
        
    def add_case(self, case, path=None, label=None): # Dictionary can be overwritten
        # Add a SatComp_case object to the cases dictionary
        if path == None:
            self.cases[case] = SatComp_case(self.case_dir,case)
        else:
            self.cases[case] = SatComp_case(path, case)
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
        obs_source, obs_label = self.__get_data_source(var)
        if not bias:
            if season:
                _season_da = season_mean(obs_source[var]).where(np.absolute(obs_source['lat'])<82)
                obs_da = _season_da[self.seas_dict[season]].to_dataset(name=var)
            self.__standard1Dplot(var, obs_da, axes, bias=False, lat_range=lat_range, 
                                  label=obs_label)
#             self.__standard1Dplot(var, obs_source, axes,
#                                   lat_range=lat_range, label=obs_label)
        
        for k in self.cases:
            _run = self.cases[k]
            _da = _run.case_da
            if season:
                _season_da = season_mean(_da[var]).where(np.absolute(_da['lat'])<82)
                _da = _season_da[self.seas_dict[season]].to_dataset(name=var)
            
            # Works with percents only
            self.__standard1Dplot(var, _da, axes, bias=bias, lat_range=lat_range, label=_run.label)
            
        fig.legend()
        
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
                                      label=obs_label)
            
            for k in self.cases:
                _run = self.cases[k]
                _da = _run.case_da
                if season:
                    _season_da = season_mean(_da[var]).where(np.absolute(_da['lat'])<82)
                    _da = _season_da[self.seas_dict[season]].to_dataset(name=var)
                
                self.__standard1Dplot(var, _da, ax, bias=bias,
                                      lat_range=lat_range,label=_run.label)
        
        fig.legend(labels=labels)
        
        xlabels = varlist; xax = axes
        if bias:
            self.add_labels(xlabels=xlabels, xaxes=xax, ylabels=["Bias (Model-Observation)"], 
                            yaxes=[axes[0]])
        else:
            self.add_labels(xlabels=xlabels, xaxes=xax)
        self.share_ylims(axes)
        
        fig.subplots_adjust(hspace=0.5)
        
        return fig
        
    def __seasonal1Dplotwrapper(self, var, bias=False, lat_range=None, **kwargs):
        '''Plot variable against latitude for observations and loaded model runs.'''
        obs_source, obs_label = self.__get_data_source(var)
        fig, axes = plt.subplots(nrows=1,ncols=4,figsize=[15,4])
    
        if bias: 
            labels = self.case_labels
        else:
            labels = [obs_label] + self.case_labels
            _season_da = season_mean(obs_source[var]).where(np.absolute(obs_source['lat'])<82)
            for season,ax in zip(self.seasons,axes):
                seas = _season_da.sel(season=season)
                self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
                                      bias=False, lat_range=lat_range)
            
#             for seas,ax in zip(_season_da,axes): # jks
#                 self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
#                                       bias=False, lat_range=lat_range)
            
        # Iterate cases first, and then seasons within that
        for k in self.cases:
            _run = self.cases[k]
            _da = _run.case_da
            
            _season_da = season_mean(_da[var])
            for season,ax in zip(self.seasons,axes):
                seas = _season_da.sel(season=season)
#                 print(seas.season.values)
                self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
                                      bias=bias,lat_range=lat_range, label=_run.label)
                
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
        obs_source, obs_label = self.__get_data_source(varlist[0])
        
        # Preprocessing step:
        for var, xax in zip(varlist, axes): # iterate over rows (i.e. layer vars)
            var_processed = []
            if bias: 
                labels = self.case_labels
            else:
                labels = [obs_label] + self.case_labels
                _season_da = season_mean(obs_source[var]).where(np.absolute(obs_source['lat'])<82)

                for season,ax in zip(self.seasons,axes):
                    seas = _season_da.sel(season=season)
                    self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
                                          bias=False, lat_range=lat_range, label=obs_label)
                
#                 for seas,ax in zip(_season_da,xax): #iterate over columns (seasons)
#                     self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
#                                           bias=False, lat_range=lat_range, label=obs_label)
            
            for k in self.cases:
                _run = self.cases[k]
                _da = _run.case_da
                _season_da = season_mean(_da[var])
                
                for season,ax in zip(self.seasons,axes):
                    seas = _season_da.sel(season=season)
                    self.__standard1Dplot(var, seas.to_dataset(name=var),ax,
                                          bias=bias, lat_range=lat_range, label=_run.label)
                
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
            val = val - obs

            im = val.sel(lat=slice(lat_lims[0],lat_lims[1])).plot(ax=ax, add_legend=False, **kwargs)
            ax.hlines(0, lat_lims[0], lat_lims[1], colors='gray', linestyles='dashed', label='')
        else:
            im = val.sel(lat=slice(lat_lims[0],lat_lims[1])).plot(ax=ax, add_legend=False, **kwargs)
            
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_title('')
        
    def plot2D(self, var, projection='PlateCarree',layers=False,
               seasonal=False, season=None, bias=False,**kwargs):
        '''
        General surface plot function. Calls more specific plotting functions.
        Supports projections:
        PlateCarree, NorthPolarStereo, Mollweide, SouthPolarStereo
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
            except:
                print('''Layers prefix %s not found. Select from: \n
                %s ''' % (var, str(list(self.layer_prefixes.keys()))))
                return None
            out = self.__layersplotwrapper(varlist, proj, bias=bias, season=season,**kwargs)
        elif seasonal:
            out = self.__seasonalplotwrapper(var, proj, bias=bias,**kwargs)
        else:
            out = self.__standard2Dplotwrapper(var, proj, bias=bias, season=season,**kwargs)
                
        return out
        
    def __standard2Dplotwrapper(self, var, projection, bias=False, season=None, **kwargs):
        '''
        Create subplots object and iterate over cases to plot with correct projection, etc.
        '''
        obs_source, obs_label = self.__get_data_source(var)
        if bias:
            fig, axes = sp_map(nrows=len(self.cases), ncols=1,
                           projection=projection, figsize=[15,2*(len(self.cases))])
            ylabels = self.case_labels
            _axes = axes
            
        else:
            fig, axes = sp_map(nrows=len(self.cases)+1, ncols=1,
                           projection=projection, figsize=[15,2*(len(self.cases)+1)])
            ylabels = [obs_label] + self.case_labels
            
            # jks handle season
            obs_da = obs_source
            if season:
                _season_da = season_mean(obs_source[var]).where(np.absolute(obs_source['lat'])<82)
                obs_da = _season_da[self.seas_dict[season]].to_dataset(name=var)
            self.standard2Dplot(var, obs_da, axes.flat[0],
                                projection=projection, bias=False)
            _axes = axes.flat[1:]
        
        for ax, k in zip(_axes, self.cases):
            _run = self.cases[k]
            _da = _run.case_da
            
            # Works with percents only
            if season:
                _season_da = season_mean(_da[var])
                _da = _season_da[self.seas_dict[season]].to_dataset(name=var)
                
            _ax, _im = self.standard2Dplot(var,_da, ax,projection=projection,
                                           bias=bias)
                        
        yax = axes
        xlabels = [var]; xax = [axes[0]]
        
        self.add_labels(ylabels=ylabels, yaxes=yax, xlabels=xlabels, xaxes=xax)
        
        self.share_clims(fig)
        
        cbar = fig.colorbar(_im, ax=axes.ravel().tolist()) # this doesn't get the right bounds
        if bias:
            cbar.set_label("Bias (Model - Observations) of %s (%s)" % (_da[var].long_name,_da[var].units), fontsize=12)
#             cbar.set_label("Bias (%s)" % _da[var].units)
        else:
            cbar.set_label("%s (%s)" % (_da[var].long_name,_da[var].units), fontsize=12)

        plt.show()
        return fig

    def __layersplotwrapper(self, varlist, projection, bias=False, season=None, **kwargs):
        '''
        Create subplots object and iterate over cases to plot with correct projection, etc.
        '''
        obs_source, obs_label = self.__get_data_source(varlist[0]) # layers must share a source
        
        if bias:
            fig, axes = sp_map(nrows=len(self.cases), ncols=len(varlist),
                           projection=projection, figsize=[2.5*len(varlist),2*(len(self.cases))])
            ylabels = self.case_labels
            
        else:
            fig, axes = sp_map(nrows=len(self.cases)+1, ncols=len(varlist),
                           projection=projection, figsize=[3*len(varlist),2*(len(self.cases)+1)])
            ylabels = [obs_label] + self.case_labels
            
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
                                projection=projection, bias=False,**kwargs)
                case_ax = varax[1:]
                
        
            for ax, k in zip(case_ax, self.cases):
                _run = self.cases[k]
                _da = _run.case_da
            
                if season:
                    _season_da = season_mean(_da[var])
                    _da = _season_da[self.seas_dict[season]].to_dataset(name=var)

                # Works with percents only
                _ax, _im = self.standard2Dplot(var,_da, ax,projection=projection,
                                               bias=bias,**kwargs)
        
        yax = axes[:,0]
        xlabels = []
        for i in varlist: # jks improve labelling with self.var_label_dict
            try:
                xlabels.append(self.var_label_dict[i])
            except:
                xlabels.append(i)
                
        xax = axes[0,:]
        self.add_labels(ylabels=ylabels, yaxes=yax, xlabels=xlabels, xaxes=xax)

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
            ylabels = self.case_labels
            _axes = axes
            
        else:
            fig, axes = sp_map(nrows=len(self.cases)+1, ncols=4,
                           projection=projection, figsize=[10,2*(len(self.cases)+1)])
            ylabels = [obs_label] + self.case_labels
            
            _season_da = season_mean(obs_source[var]).where(np.absolute(obs_source['lat'])<82)

            for season,ax in zip(self.seasons,axes[0]):
                seas = _season_da.sel(season=season)
                _ax, _im = self.standard2Dplot(var, seas.to_dataset(name=var), ax,
                                    projection=projection, bias=False)
                
#             for seas,ax in zip(_season_da,axes[0]):
#                 _ax, _im = self.standard2Dplot(var, seas.to_dataset(name=var), ax,
#                                     projection=projection, bias=False)
            _axes = axes[1:]
            
        # Iterate cases first, and then seasons within that
        for xaxes, k in zip(_axes, self.cases):
            _run = self.cases[k]
            _da = _run.case_da
            
            _season_da = season_mean(_da[var])
            
            for season,ax in zip(self.seasons,xaxes):
                seas = _season_da.sel(season=season)
                _ax, _im = self.standard2Dplot(var, seas.to_dataset(name=var), ax,
                                projection=projection, bias=bias)
            
#             for seas,ax in zip(_season_da,xaxes):
#                 _ax, _im = self.standard2Dplot(var, seas.to_dataset(name=var), ax,
#                                 projection=projection, bias=bias)

        self.share_clims(fig)
        
        cbar = fig.colorbar(_im, ax=axes.ravel().tolist())
        if bias:
            cbar.set_label("Bias (%s)" % _da[var].units)
        else:
            cbar.set_label("%s (%s)" % (_da[var].long_name,_da[var].units), fontsize=12)
        
        yax = axes[:,0]
        xlabels = self.seasons; xax = axes[0,:]
        self.add_labels(ylabels=ylabels, yaxes=yax, xlabels=xlabels, xaxes=xax)
                
        return fig
    
    def standard2Dplot(self, var, da, ax, projection=None, bias=False, **kwargs):
        '''
        Basic workhorse plotting function. needs to handle bias and seasons
        '''
        obs_source, obs_label = self.__get_data_source(var)
        lat_lims = [-90,90]
        if projection == ccrs.NorthPolarStereo(): 
            lat_lims = [59.5,90]
            polarCentral_set_latlim(lat_lims, ax)
        if projection == ccrs.SouthPolarStereo(): 
            lat_lims = [-90,-60]
            polarCentral_set_latlim(lat_lims, ax)
        if 'time' not in da.dims: # catches a bug with seasonal averages
            val = da[var].where(da['lat'] > lat_lims[0])
        else:
            val = da[var].mean(dim = 'time', skipna=True).where(
                                da['lat'] > lat_lims[0])            
        
        if bias:
            try: # get season observations if passed da is seasonally processed
                season = val['season'].values
                obs = season_mean(obs_source[var]).sel(season=season).where(
                                  (da['lat'] > lat_lims[0]) & (np.absolute(obs_source['lat'])<82))

            except:
                obs = obs_source[var].mean(dim = 'time', skipna=True).where(
                                                da['lat'] > lat_lims[0])
            obs = obs.interp_like(val) # quick interp fix for weird grid mismatch (bad.)
#             val = val.interp_like(obs) # quick interp fix for weird grid mismatch (bad.)

            val = val - obs
    
    

#             im = val.plot(ax=ax,cmap=plt.get_cmap('bwr'),transform=ccrs.PlateCarree(),
#                           add_colorbar=False, robust=True, **kwargs) #jks robust good?
#             im = val.plot(ax=ax,cmap=plt.get_cmap('icefire'),transform=ccrs.PlateCarree(),
#                           add_colorbar=False, robust=True, **kwargs) #jks robust good?
            if 'cmap' in kwargs.keys():
                im = val.plot(ax=ax,transform=ccrs.PlateCarree(),
                              add_colorbar=False, robust=True, **kwargs)
            else:
                im = val.plot(ax=ax,cmap=plt.get_cmap('bwr'),transform=ccrs.PlateCarree(),
                              add_colorbar=False, robust=True, **kwargs)
                
        else:
            im = val.plot(ax=ax,cmap=plt.get_cmap('jet'),transform=ccrs.PlateCarree(),
                      add_colorbar=False, robust=True, **kwargs) # robust good?
            
        ax.set_title('')
#             im = val.plot(ax=ax,cmap=plt.get_cmap('jet'),transform=ccrs.PlateCarree(),
#                       add_colorbar=False, vmin=0, vmax=100, **kwargs)


        add_map_features(ax)
                    
        return ax, im
    
    def band_biases(self, var):
        '''
        Calculate the model bias with respect to observations over 10 degree(ish) latitude bands.
        '''
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
        fig = plt.figure(figsize=[6*(len(self.cases)),10],dpi=200)
#         months = ['J','F','M','A','M','J','J', 'A','S','O','N','D'] # month initials
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
        
    def add_labels(self, xlabels=None, ylabels=None, xaxes=None, yaxes=None):
        try:
            for ax,label in zip(xaxes,xlabels):
                ax.text(0.5, 1.1, label, va='bottom', ha='center',
                    rotation='horizontal', rotation_mode='anchor',
                    transform=ax.transAxes, fontsize=14)
        except: pass

        try:
            for ax,label in zip(yaxes,ylabels):
                ax.text(-0.2, 0.5, label, va='bottom', ha='center',
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
                if isinstance(j, mpl.collections.QuadMesh):
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
        if var in self.goccp_data.data_vars:
            data_source = self.goccp_data
            label = "GOCCP"
        elif var in self.ceres_data.data_vars:
            data_source = self.ceres_data
            label = "CERES-EBAF"
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
        
class SatComp_case:
    
    def __init__(self, casedir, case):
        self.case_dir = casedir # General directory where all cases are stored
        self.case = case # The case name
        
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
        
    def add_ds(self):
        try:
            print('Trying to load concatenated file for %s' % self.case)
            _ds = xr.open_dataset('%s%s/%s.nc' % (self.case_dir, self.case, self.case))
        except: # this probably wouldn't work anyway, at least on NIRD
            print('Failed, using xr.open_mfdataset.')
            self.geth0s() # create list of appropriate output files self.files
            _ds = xr.open_mfdataset(self.files)#, combine='by_coords') #, chunks={'lat':10}) #chunk?
        # Shift time-bounds to fix comparison with satellites (from Jen Kay)
        _ds['time'] = _ds['time_bnds'].isel(bnds=0)
            
        try:
            ds = _ds.sel(time=slice('2009-06-01', '2013-05-01')) # gets all 48 files
        except:
            print("Four-year format not found.")
        ds = add_weights(ds)
        
        # Normalize calculated SLF values.
        ds['IC_SLF'] = ds['SLFXCLD_ISOTM']/ds['CLD_ISOTM']
        ds['CT_SLF'] = ds['CT_SLFXCLD_ISOTM']/ds['CT_CLD_ISOTM']

        self.case_da = ds
        
        print("%s load successfully." % self.case)
        
    def geth0s(self):
        prepath = '%s%s/atm/hist/' % (self.case_dir, self.case)
        allfiles = os.listdir(prepath)
        self.files = [prepath + i for i in allfiles if "h0" in i]
    