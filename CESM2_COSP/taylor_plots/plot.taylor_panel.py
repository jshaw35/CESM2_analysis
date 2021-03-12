#!/usr/bin/env python

# Original Code by Hillman
# Some comments and changes from Jonah Shaw

import numpy as numpy
import Scientific.IO.NetCDF as NetCDF
import taylor as taylor
import matplotlib.pyplot as pyplot
import matplotlib as matplotlib
import matplotlib.patches as patches

def get_data(inptr,varname):
    '''returns an array where the _FillValues or 
    missing_value cells are masked, great!'''
    
    if '_FillValue' in dir(inptr.variables[varname]):
        data = numpy.ma.masked_values(
                inptr.variables[varname].getValue(),
                getattr(inptr.variables[varname],'_FillValue'),
            )
    elif 'missing_value' in dir(inptr.variables[varname]):
        data = numpy.ma.masked_values(
                inptr.variables[varname].getValue(),
                getattr(inptr.variables[varname],'missing_value'),
            )
    else:
        data = inptr.variables[varname].getValue()
    return data

# Control names dictionary
cntlnames = {
        'SWCFTOA': 'CERES-EBAF',
        'LWCFTOA': 'CERES-EBAF',
        'CLDTOT_ISCCPCOSP': 'ISCCP',
        'CLDTOT_MISR': 'MISR',
        'CLDTOT_CAL': 'CALIPSO',
        'CLDLOW_CAL': 'CALIPSO',
        'CLDMED_CAL': 'CALIPSO',
        'CLDHGH_CAL': 'CALIPSO',
        'CLDLOW_THICK_MISR': 'MISR',
        'CLDHGH_THICK_MODIS': 'MODIS',
    }

# Case names
testnames = ('CAM4','CAM5')
testcolors = ('SkyBlue','Firebrick')

# Cloud forcing plot
varlist = (
        ('SWCFTOA','LWCFTOA'),
        ('CLDTOT_ISCCPCOSP','CLDTOT_MISR','CLDTOT_CAL'),
        ('CLDLOW_CAL','CLDMED_CAL','CLDHGH_CAL'),
        ('CLDLOW_THICK_MISR','CLDHGH_THICK_MODIS'),
    )

# Fix fonts
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['text.usetex'] = True

# Open figure
figure = pyplot.figure()

for iplot,varnames in enumerate(varlist):
    print varnames
    nvars = len(varnames)
    ntest = len(testnames)

    cc = numpy.zeros([nvars,ntest])
    ratio = numpy.zeros([nvars,ntest])
    bias = numpy.zeros([nvars,ntest])

    # Loop over variables
    for ivar,varname in enumerate(varnames):

        # read data, from the control file (obs?)
        cntl_file = 'plotvars/'+cntlnames[varname]+'.'+varname+'.nc'
        cntl_inptr = NetCDF.NetCDFFile(cntl_file,'r')
        cntl_data = get_data(cntl_inptr,varname) # requested variable field
        cntl_lat = cntl_inptr.variables['lat'].getValue() # lat coord
        cntl_lon = cntl_inptr.variables['lon'].getValue() # lon coord

        # Loop over test cases
        for itest,testname in enumerate(testnames):

            print '    variable: '+varname+', case: '+testname

            # Read data, from model
            test_file = 'plotvars/'+testname+'.'+varname+'.nc'
            test_inptr = NetCDF.NetCDFFile(test_file,'r')
            test_data = get_data(test_inptr,varname)
            test_lat = test_inptr.variables['lat'].getValue()
            test_lon = test_inptr.variables['lon'].getValue()

            # Calculate and save statistics, JKS nice this is an easy call to make.
            stats = taylor.Taylor_statistics(
                    test_data,test_lat,test_lon,
                    cntl_data,cntl_lat,cntl_lon,
                )
            cc[ivar,itest] = stats.cc
            ratio[ivar,itest] = stats.ratio
            bias[ivar,itest] = stats.bias

            #print 'cc = %.2f, ratio = %.2f, bias = %.2f'%(stats.cc,stats.ratio,stats.bias)

    # Make plot
    ax = figure.add_subplot(2,2,iplot+1,frameon=False)
    taylor_diagram = taylor.Taylor_diagram(
            ax,cc,ratio,bias,
            casecolors=testcolors,
            varlabels=range(1,nvars+1),
        )

    # Reference bias bubbles, wut is this?
    ref_bias = 0.1 # This is a 10% bias reference bubble in the lower-left corner
    yloc = 0.05*taylor_diagram.xymax + ref_bias/2.0
    xloc = 0.05*taylor_diagram.xymax + ref_bias/2.0
    circle = patches.Circle(
            (xloc,yloc),ref_bias/2.0,
            color="black",
            alpha=0.30,
        )
    ax.add_patch(circle)

    # Reference bias bubble points - centered at the reference bubble
    circle = patches.Circle(
            (xloc,yloc),0.01,
            color="black",
        )
    ax.add_patch(circle)

    # Reference bias text
    ax.text(
            xloc+ref_bias/2.0 + 0.01*taylor_diagram.xymax,yloc,
            "%.0f%s bias"%(ref_bias*100,r"\%"),
            color="Black",
            fontsize=8,
            horizontalalignment="left",
            verticalalignment="center"
        )

    # Case labels
    xloc = taylor_diagram.xymax*0.95
    yloc = taylor_diagram.xymax*0.05
    dy = taylor_diagram.xymax*0.05
    for itest,testname in enumerate(testnames[::-1]):
        ax.text(
                xloc,yloc+itest*dy, # place these just above the dots
                testname,
                color=testcolors[::-1][itest],
                fontsize=8,
                horizontalalignment="right",
                verticalalignment="bottom",
            )

figure.savefig('plots/taylor_panel.pdf',format='pdf')
