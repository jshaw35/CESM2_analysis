from Common.imports import *
from functions import *

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy as cy

def org_feedbacks_case(case_path):
    '''
    This function take a directory path that 
    contains feedback files (via Lise script).
    It opens each feedback file and extracts the needed feedback variables,
    which are stored in a single dataset object (returned).
    '''
    fb_dir ='feedbacks_soden/mod_by_me/output_jonahks/'
    _fbs = os.listdir("%s/%s" % (fb_dir,case_path)) # List feedback files
    
    # Calculate the averaged Arctic quantities for each file appropriately:
    # Albedo feedback
    _albedo_file = [i for i in _fbs if "albedo_feedback" in i][0]
    _albedo = add_weights(xr.open_dataset("%s/%s/%s" % (fb_dir,case_path,_albedo_file)).sel(run=0))
    _albedo_vars = ['fb_albedo_ttsky']
    
    # Create a new dataset object to store specific variables within
    _all_fb = _albedo[_albedo_vars[0]].to_dataset(name=_albedo_vars[0])
    _albedo.close()
    
    # Cloud feedbacks
    _cld_file = [i for i in _fbs if "cloud_feedback" in i][0] # select file, a little clunky to avoid using regular expressions
    _cld = add_weights(xr.open_dataset("%s/%s/%s" % (fb_dir,case_path,_cld_file)).sel(run=0))    
    cld_vars = ["fb_cloud_sw","fb_cloud_lw","fb_adj_cloud_sw","fb_adj_cloud_lw"]
    for _var in cld_vars:
        _all_fb[_var] = _cld[_var]
    _cld.close()
    
    # Water Vapor feedback (separate files)
    _wv_sw_file = [i for i in _fbs if "water_vapor_sw" in i][0]    
    _wv_lw_file = [i for i in _fbs if "water_vapor_lw" in i][0]
    _wv_sw = add_weights(xr.open_dataset("%s/%s/%s" % (fb_dir,case_path,_wv_sw_file)).sel(run=0))
    _wv_lw = add_weights(xr.open_dataset("%s/%s/%s" % (fb_dir,case_path,_wv_lw_file)).sel(run=0))
    _wv_vars = ["fb_water_vapor_sw_ttsky", "fb_water_vapor_lw_ttsky"] # Will need to zip
    for _var,_da in zip(_wv_vars,[_wv_sw,_wv_lw]):
        _all_fb[_var] = _da[_var]
        _da.close()
        
    # Lapse Rate feedback
    _lr_file = [i for i in _fbs if "lapse_rate" in i][0]
    _lr = add_weights(xr.open_dataset("%s/%s/%s" % (fb_dir,case_path,_lr_file)).sel(run=0))
    _lr_vars =["fb_lapse_rate_ttsky"]
    for _var in _lr_vars:
        _all_fb[_var] = _lr[_var]
    _lr.close()
    
    # Planck feedback
    _planck_file = [i for i in _fbs if "planck_feedback" in i][0]
    _planck = add_weights(xr.open_dataset("%s/%s/%s" % (fb_dir,case_path,_planck_file)).sel(run=0))
    _planck_vars = ['fb_planck_ttsky']
    for _var in _planck_vars:
        _all_fb[_var] = _planck[_var]
    _planck.close()

    _all_fb = add_weights(_all_fb)
    
    return _all_fb

def comp_feedbacks(fb_ds):
    '''
    This function averages feedbacks in the Arctic Circle,
    weighting by the gridcell area and month duration.
    '''
    # Used for weighting by month. I tested this and it didn't make a big change here, but could elsewhere.
    wgt_mon=[31,28,31,30,31,30,31,31,30,31,30,31]
    month_length = xr.DataArray(wgt_mon, coords=[fb_ds['fb_albedo_ttsky']['month']], name='month_length')

    # This matrix cross-product will allow me to do all of my weighting in a single step.
    all_weights = month_length @ fb_ds['cell_weight']
    
    for _var in fb_ds.data_vars:
        _ds = fb_ds[_var]
        mask = _ds['lat'] < 66
        _arc_val = masked_average(_ds,dim=['lat','lon','month'],weights=all_weights,mask=mask)
        print(_var, ": ", _arc_val.values)
        

def autolabel2(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('%.2f' %height,
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    fontsize=11,
                    ha='center', va='bottom')
        

def autolabel3(rects,ax,case=None):
    """Attach a text label above each bar in *rects*, displaying its height.
    Case is just an exception for plotting negative values on a flipped axis. 
    """
    for rect in rects:
        height = rect.get_height()
        if (height < 0 and not case): # for consistent placement
            txtxy = (0,-14)
        else:
            txtxy = (0,3)
        ax.annotate('%.2f' %height,
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=txtxy,  # 3 points vertical offset
                    textcoords="offset points",
                    fontsize=11,
                    ha='center', va='bottom')
        
def autolabel4(line,ax,case=None):
    """
    Autolabel but for normal plots. In development. 
    """
    for _x,_y in zip(line.get_xdata(),line.get_ydata()):
        if (_y < 0 and not case): # for consistent placement
            txtxy = (0,-16)
        else:
            txtxy = (0,5)
        ax.annotate('%.2f' %_y,
                    xy=(_x, _y),
#                     xy=(rect.get_x() + rect.get_width() / 2, _y),
                    xytext=txtxy,  # 3 points vertical offset
                    textcoords="offset points",
                    fontsize=11,
                    ha='center', va='bottom')
        
        
def weight_feedbacks(case_dict,weights=None,labels=None):
    '''
    Takes a dictionary of organized feedbacks/cases.
    Computes an annual average for each feedback and
    stores in individual lists in a returned dictionary object.
    Currently adding functionality to handle unnormalized feedbacks.
    Weights will be a dictionary with identical keys to case_dict.
    '''
    
    # Used for weighting by month. I tested this and it didn't make a big change here, but could elsewhere.
    wgt_mon=[31,28,31,30,31,30,31,31,30,31,30,31]
    
    # Pick a random variable to use coords from here
    sample = next(iter( case_dict.values() ))
#     sample = case_dict[case_dict.keys()[0]]['fb_albedo_ttsky']['month']
    month_length = xr.DataArray(wgt_mon, coords=[sample['fb_albedo_ttsky']['month']], name='month_length')

    # This matrix cross-product will allow me to do all of my weighting in a single step.
    all_weights = month_length @ sample['cell_weight']
    
    barvals_PL = []
    barvals_LR = []
    barvals_WV = []
    barvals_A = []
    barvals_CLSW = []
    barvals_CLLW = []
    
    cases = []
    if labels:
        keys = labels.keys()
    else: 
        keys = case_dict.keys()
    for _mod in keys:
        print(_mod)
#         print(weights[_mod])
        _da = case_dict[_mod]

        _avg_dict = {}
        cases.append(_mod)
        for _var in _da.data_vars:
            _ds = _da[_var]
            mask = _ds['lat'] < 66
            _arc_val = masked_average(_ds,dim=['lat','lon','month'],weights=all_weights,mask=mask)
    #             print(_var, ": ", _arc_val.values)
            if weights: # Normalize by Arctic increase
                _avg_dict[_var] = _arc_val.values / weights[_mod]
            else:
                _avg_dict[_var] = _arc_val.values 

        barvals_PL.append(_avg_dict['fb_planck_ttsky'])        
        barvals_LR.append(_avg_dict['fb_lapse_rate_ttsky'])
        barvals_WV.append(_avg_dict['fb_water_vapor_sw_ttsky']+_avg_dict['fb_water_vapor_lw_ttsky'])
        barvals_A.append(_avg_dict['fb_albedo_ttsky'])
        barvals_CLSW.append(_avg_dict['fb_adj_cloud_sw'])
        barvals_CLLW.append(_avg_dict['fb_adj_cloud_lw'])
            
    outdict = dict(planck_fb = barvals_PL,
                   lapserate_fb = barvals_LR,
                   watervapor_fb = barvals_WV,
                   albedo_fb = barvals_A,
                   cloudsw_fb = barvals_CLSW,
                   cloudlw_fb = barvals_CLLW
                  )
    
    return outdict


def barplot_feedbacks1(case_dict,weighted_fbs,labels=None):
    '''
    Plot the weighted annual feedbacks using a bar plot. All feedbacks.
    '''
    
    fig, ax = plt.subplots(figsize=(15.,6))

    if labels:
        keys = labels.values()
    else:
        keys = case_dict.keys()
    xind = range(len(keys))
    width = 0.167         # the width of the bars
    xind_p1 = [i-(width*2.5) for i in xind]
    xind_p2 = [i-(width*1.5) for i in xind]
    xind_p3 = [i-(width*0.5) for i in xind]
    xind_p4 = [i+(width*0.5) for i in xind]
    xind_p5 = [i+(width*1.5) for i in xind]
    xind_p6 = [i+(width*2.5) for i in xind]
    p1 = ax.bar(xind_p1, weighted_fbs['planck_fb'], width, label='Planck')
    p2 = ax.bar(xind_p2, weighted_fbs['lapserate_fb'], width, label='Lapse Rate')
    p3 = ax.bar(xind_p3, weighted_fbs['watervapor_fb'], width, label='Water Vapor')
    p4 = ax.bar(xind_p4, weighted_fbs['albedo_fb'], width, label='Albedo')
    p5 = ax.bar(xind_p5, weighted_fbs['cloudsw_fb'], width, label='adj. SW-cloud')
    p6 = ax.bar(xind_p6, weighted_fbs['cloudlw_fb'], width, label='adj. LW-cloud')

    autolabel3(p1,ax)
    autolabel3(p2,ax)
    autolabel3(p3,ax)
    autolabel3(p4,ax)
    autolabel3(p5,ax)
    autolabel3(p6,ax)

    ax.set_xticks(xind)
    if labels:
        ax.set_xticklabels(keys, fontsize=12,rotation=30)
    else:
        ax.set_xticklabels(keys, rotation=45)#, fontsize=15)

    ax.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]), ('Planck', 'Lapse Rate', 'Water Vapor', 'Albedo', 'adj. SW-cloud', 'adj. LW-cloud'))
    ax.set_ylabel('[Wm$^{-2}$K$^{-1}$]',fontsize=18)
    
    return fig,ax


def barplot_feedbacks2(case_dict,weighted_fbs,labels=None):
    '''
    Plot the weighted annual feedbacks using a bar plot. Some feedbacks and the sum.
    '''
    
    fig, ax = plt.subplots(figsize=(15.,6))

    if labels:
        keys = labels.values()
    else:
        keys = case_dict.keys()
    xind = range(len(keys))

    width = 0.25         # the width of the bars
    xind_p1 = [i-(width*1.5) for i in xind]
    xind_p2 = [i-(width*0.5) for i in xind]
    xind_p3 = [i+(width*0.5) for i in xind]
    xind_p4 = [i+(width*1.5) for i in xind]
    
    fbs = [weighted_fbs['planck_fb'], weighted_fbs['lapserate_fb'],
           weighted_fbs['watervapor_fb'], weighted_fbs['albedo_fb'],
           weighted_fbs['cloudsw_fb'], weighted_fbs['cloudlw_fb']]
    
    p1 = ax.bar(xind_p1, weighted_fbs['cloudsw_fb'], width, label='adj. SW-cloud')
    p2 = ax.bar(xind_p2, weighted_fbs['cloudlw_fb'], width, label='adj. LW-cloud')
    p3 = ax.bar(xind_p3, weighted_fbs['lapserate_fb'], width, label='Lapse Rate')
    net_fb = np.sum(fbs, axis=0)
    p4 = ax.bar(xind_p4, net_fb, width, label='Net Feedback')

    autolabel3(p1,ax)
    autolabel3(p2,ax)
    autolabel3(p3,ax)
    autolabel3(p4,ax)
#     plt.ylim(min(barvals_PL)-0.05,abs(min(barvals_PL))+0.05)

    ax.set_xticks(xind)
    if labels:
        ax.set_xticklabels(keys, fontsize=12,rotation=30)
    else:
        ax.set_xticklabels(keys, rotation=45)#, fontsize=15)

    ax.legend((p1[0], p2[0], p3[0], p4[0]), ('adj. SW-cloud', 'adj. LW-cloud','Lapse Rate','Net'))

    ax.set_ylabel('[Wm$^{-2}$K$^{-1}$]',fontsize=18)
    
    return fig,ax


def calc_arc_dT(initial_ts,final_ts,lat_range=[66,90]):
    '''
    This function takes the path to inital and final 
    surface temperature files. It averages over months and
    calculates the average annual warming in the Arctic (>66N), 
    weighting by cell area and month length.
    '''
    
    # Open TS timeseriesinput files
    ts_i = xr.open_dataset(initial_ts)
    ts_f = xr.open_dataset(final_ts)

    # Take the difference and average across months.
    _d_ts = (ts_f['TS'] - ts_i['TS']).groupby('time.month').mean('time')
    d_ts = add_weights(_d_ts.to_dataset())

    # Create month length object for weighting
    wgt_mon=[31,28,31,30,31,30,31,31,30,31,30,31] # month length (days)
    month_length = xr.DataArray(wgt_mon, coords=[d_ts['month']], name='month_length')

    # Combine weights via matrix product
    all_weights = month_length @ d_ts['cell_weight']

#     mask = d_ts['lat'] < 66 # Mask below the Arctic
    mask = np.bitwise_or(d_ts['lat']<=lat_range[0],d_ts['lat']>lat_range[1])
    _arc_val = masked_average(d_ts['TS'],dim=['lat','lon','month'],weights=all_weights,mask=mask)
    
    return _arc_val
    
    
def plot_months_line(dict_cases, var, lat_range=[66,82], dTS=None,
                     bias=False, ax=None, labels=None, **kwargs):
    colors = sns.color_palette("colorblind")
    months = ['J','F','M','A','M','J','J','A','S','O','N','D'] # month initials
    
    if ax:
        axes = ax
        plt.sca(ax)
        fig=None
    else:
        fig, axes = plt.subplots(nrows=1,ncols=1,figsize=[15,10])
        
    lines = []
    for k,color in zip(dict_cases,colors):
        _run = dict_cases[k]
        mon_vals = _run[var].sel(lat=slice(lat_range[0],lat_range[1]))

        mon_vals2 = add_weights(mon_vals)
        _weights = mon_vals['cell_weight']
            
        out_dat = masked_average(mon_vals2, dim=['lat','lon'], weights=_weights)
        if dTS:
            out_dat = out_dat / dTS[k]
        if labels:
            _lbl = labels[k]
        else:
            _lbl = k

        _ln = out_dat.plot(ax=axes,label=_lbl, color=color, **kwargs)

    plt.xticks(np.arange(0,len(months)+1,1), months) # this works I guess
    plt.legend(labels.keys()) # JKS
    
    axes.set_xlabel('Month',fontsize=16)
    axes.set_title('')    
    axes.set_ylabel("%s (%s)" % (dict_cases[k][var].long_name,'Wm$^{-2}$K$^{-1}$'), fontsize=16)
    axes.hlines(0,0,12, linestyle='dashed',color='gray')
    axes.set_xbound(0,11)    

    return fig, axes


def barplot_single(case_dict,weighted_fbs,var,labels=None,vals=None, ax=None,**kwargs):
    '''
    Plot the weighted annual feedbacks across the different models 
    for a single variable.
    '''
    
    if ax:
        axes = ax
        plt.sca(ax)
        fig=None
    else:
        fig, ax = plt.subplots(figsize=(8.,6))
#         fig, axes = plt.subplots(nrows=1,ncols=1,figsize=[15,10])
    
    colors = sns.color_palette('colorblind')    
    
    if labels:
        keys = labels.values()
    else:
        keys = case_dict.keys()
    
    if vals:
        sorted_vals = []
        val_dict = {}
        for i in case_dict.keys():
            sorted_vals.append(vals[i])
            val_dict[i] = vals[i]
        p1 = ax.bar(case_dict.keys(), sorted_vals,color=colors,**kwargs) #label='adj. SW-cloud'
        val_dict[i] = vals[i]
    else:
        p1 = ax.bar(keys, weighted_fbs[var],color=colors,**kwargs)
        val_dict = dict(zip(case_dict.keys(),weighted_fbs[var]))

    autolabel3(p1,ax)

    if labels:
        ax.set_xticklabels(keys, fontsize=12,rotation=30)
    else:
        ax.set_xticklabels(keys, rotation=45)#, fontsize=15)

#     ax.legend((p1[0], p2[0], p3[0], p4[0]), ('adj. SW-cloud', 'adj. LW-cloud','Lapse Rate','Net'))
    ax.set_ylabel('[Wm$^{-2}$K$^{-1}$]',fontsize=18)
    
    return fig,ax,val_dict