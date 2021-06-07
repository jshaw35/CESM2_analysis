# Original code by Hillman.
# Updating by Jonah Shaw for CESM2/CAM6 COSP comparison.

import numpy as np
import matplotlib.artist as artist
import matplotlib.patches as patches
import matplotlib.lines as lines
import matplotlib.pylab as pylab
from matplotlib.axes import Axes

def regrid(data1,lat1,lon1,lat2,lon2):
    print('nearest neighbor interpolation starting...', end = '')
    # check inputs
    if ((len(lat1) < len(lat2)) or (len(lon1) < len(lon2))):
        print('Target grid is higher resolution than input grid.')
        raise

    shape = (len(lat2),len(lon2))
    data2 = np.ma.zeros(shape,dtype=data1.dtype,) # use same resolution as data1
#     print('data2, shape: ',data2.shape)
    if 'fill_value' in dir(data1): # dir?
        data2.fill_value = data1.fill_value # mask where the target grid is masked

    # regrid data
    # This a nearest neighbor approach. Is that ok?    
    for idx2,lon in enumerate(lon2):
        for idy2,lat in enumerate(lat2):
            idx1 = np.abs(lon1 - lon).argmin()
            idy1 = np.abs(lat1 - lat).argmin() 
#             print('idx1: ',idx1,idx1.shape)
#             print('idy1: ',idy1,idy1.shape)
            data2[...,idy2,idx2] = data1[...,idy1,idx1]
    
    print('done.')
    return data2


class Taylor_statistics():

    cc = 0.0
    ratio = 0.0
    bias = 0.0
    rmsnorm = 0.0

    def __init__(
            self,
            test_data,test_lat,test_lon,
            cntl_data,cntl_lat,cntl_lon,
        ):

        # regrid data
        if len(test_lat) < len(cntl_lat):
            lat = test_lat
        else:
            lat = cntl_lat
        if len(test_lon) < len(cntl_lon):
            lon = test_lon
        else:
            lon = cntl_lon
            
        
        test_data = regrid(test_data,test_lat,test_lon,lat,lon)
        cntl_data = regrid(cntl_data,cntl_lat,cntl_lon,lat,lon)

        # calculate weights
        lat_weights = np.cos(lat[:]*np.pi/180.)
        lat_weights = lat_weights/lat_weights.sum()

        # make full lat,lon dimensioned weight array
        weight_array = np.ma.ones([len(lat),len(lon)])
        for idx in range(len(lon)):
            weight_array[:,idx] = lat_weights[:]

        # mask missing values
        test_data.mask = np.ma.mask_or(test_data.mask,cntl_data.mask)
        cntl_data.mask = np.ma.mask_or(test_data.mask,cntl_data.mask)
        weight_array.mask = test_data.mask

        # calculate statistics
        self.calculate(test_data,cntl_data,weight_array)

    def calculate(self,test,cntl,wgt):
        """Calculate Taylor statistics for making taylor diagrams."""

        # calculate sums and means
        sumwgt = np.ma.sum(wgt)
        meantest = np.ma.sum(wgt*test)/sumwgt
        meancntl = np.ma.sum(wgt*cntl)/sumwgt

        # calculate variances
        stdtest = (np.ma.sum(wgt*(test-meantest)**2.0)/sumwgt)**0.5
        stdcntl = (np.ma.sum(wgt*(cntl-meancntl)**2.0)/sumwgt)**0.5

        # calculate correlation coefficient
        ccnum = np.ma.sum(wgt*(test-meantest)*(cntl-meancntl))
        ccdem = sumwgt*stdtest*stdcntl
        self.cc = ccnum/ccdem

        # calculate variance ratio
        self.ratio = stdtest/stdcntl

        # calculate bias
        self.bias = (meantest - meancntl)/np.abs(meancntl)
        #self.bias = meantest - meancntl

        # calculate centered pattern RMS difference
        rmssum = np.ma.sum(wgt*((test-meantest)-(cntl-meancntl))**2.0)
        rmserr = (rmssum/sumwgt)**0.5
        self.rmsnorm = rmserr/stdcntl


class Taylor_diagram():
    """
    Custom axes class to facilitate making Taylor diagrams. Adds the taylor
    method to draw the diagram.
    """

    def __init__(
            self,ax,cc,ratio,bias,
            rms=None,casecolors=None,varlabels=None,
        ):

        # default casecolors
        if casecolors is None:
            casecolors = [
                    'blue','red','green','cyan','magenta','yellow',
                    '0.75','0.50','0.25','0.10',
                ]

        # copy axes instance
        self.ax = ax

        # plot size
        #self.xymax = 1.50
#         self.xymax = np.max([1.50,np.max(ratio)+np.max(np.abs(bias))/2.0]) # JKS added the np.abs() to fix errors
        self.xymax = np.max([1.50,np.max(ratio+np.abs(bias)/2.0)]) # JKS added the np.abs() and moved np.max outward to fix errors
#         print('self.xymax: ',self.xymax)
        
        self.xymax = np.ceil(self.xymax*4)/4 # JKS Round up to the nearest 0.25

#         print('self.xymax: ',self.xymax)

        # draw axes
        self.draw_axes(ratio,cc)

        # add some reference arcs
        self.draw_ratio(1.0,linestyle='solid',color='black')
        self.draw_ratio(np.min(ratio),linestyle='dotted',color='black')
        self.draw_ratio(np.max(ratio),linestyle='dotted',color='black')

        # draw some reference lines
        self.draw_radii(np.min(cc))
        self.draw_radii(np.max(cc))

        # draw points
        self.draw_point(
                ratio,cc,bias,
                bubblecolors=casecolors,
                varlabels=varlabels,
            )

        # draw rms circles
        if rms is not None:
            self.draw_rms(min(rms))
            self.draw_rms(max(rms))

    def polar_transform(self,r,theta):
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        return x,y

    def draw_radii(self,cc,color='black',linestyle='dotted'):
        theta = np.arccos(cc)
        x,y = self.polar_transform(self.xymax,theta)
        self.ax.plot([0,x],[0,y],':k',linewidth=0.5)

    def draw_ratio(self,r,color='black',linestyle='dotted'):
        arc = patches.Arc(
                [0,0],2*r,2*r,
                theta1=0.0,
                theta2=90.0,
                color=color,
                linestyle=linestyle,
                linewidth=0.5,
            )
        self.ax.add_patch(arc)

    def draw_axes(self,ratio,cc):

        # draw axes
        wedge = patches.Wedge(
                [0,0],self.xymax,
                theta1=0.0,
                theta2=90.0,
                color='black',
                fill=False
            )
        self.ax.add_patch(wedge)

        # correlation axis label
        self.ax.text(
                self.xymax*0.8,self.xymax*0.8,'Correlation',
                color='black',
                rotation='-45',
                horizontalalignment='center',
                verticalalignment='center'
            )

        # x-axis label
#         self.ax.xaxis.set_label_text('Relative standard deviation')
        # y-axis label
        self.ax.yaxis.set_label_text('Normalized standard deviation')

        # x-axis tickmarks
        xmajorticks = np.arange(0.0,self.xymax+0.01,0.25)
#         xmajorticks = np.arange(0.0,self.xymax+0.25,0.25) # JKS fix mismatch with self.xymax by rounding up in 0.25 increments
        xmajorticklabels = []
        xminorticks = [np.min(ratio),np.max(ratio)]
        xminorticklabels = ["%.2f"%(np.min(ratio)),"%.2f"%(np.max(ratio))]

        # y-axis tickmarks
        ymajorticks = np.arange(0.0,self.xymax+0.01,0.25)
#         ymajorticks = np.arange(0.0,self.xymax+0.25,0.25) # JKS fix mismatch
        ymajorticklabels = []
        for idy in range(len(ymajorticks)):
            if idy%2 == 0:
                ymajorticklabels.append('%.1f'%(ymajorticks[idy]))
            else:
                ymajorticklabels.append('')

        # correlation tickmarks
        ccmajorticks = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        ccminorticks = [0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99]
        ccmajorticklabellocs = [
                0.5,0.95,np.min(cc),np.max(cc)
            ]
        ccmajorticklabels = [
                "0.5","0.95","%.2f"%(np.min(cc)),"%.2f"%(np.max(cc))
            ]

        # Set axis limits
        self.ax.set_ylim([min(ymajorticks)-0.01, max(ymajorticks)+0.01])
#         self.ax.set_ylim([min(ymajorticks)-0.01, max(ymajorticks)*1.05]) # JKS adjustment here
        self.ax.set_xlim([min(xmajorticks)-0.01, max(xmajorticks)+0.01])
#         self.ax.set_xlim([min(xmajorticks)-0.01, max(xmajorticks)*1.05])
#         print('ymajorticks: ',ymajorticks)
#         print('xmajorticks: ',xmajorticks)
        
        # set and draw tickmarks
        self.ax.set_aspect('equal')
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.yaxis.set_ticks_position('left')
#         self.ax

        self.ax.xaxis.set_ticks(xmajorticks)
        self.ax.xaxis.set_ticklabels(xmajorticklabels)
        self.ax.xaxis.set_ticks(xminorticks,minor=True)
        self.ax.xaxis.set_ticklabels(xminorticklabels,minor=True)

        self.ax.yaxis.set_ticks(ymajorticks)
        self.ax.yaxis.set_ticklabels(ymajorticklabels)


        # draw correlation tickmarks
        majorticklength = self.xymax*0.02
        minorticklength = self.xymax*0.01
        for cctick in ccmajorticks:
            r = [self.xymax-majorticklength,self.xymax            ]
            theta = [np.arccos(cctick),np.arccos(cctick)]
            x,y = self.polar_transform(r,theta)
            self.ax.plot(x,y,'-k')
        for i,ccticklabel in enumerate(ccmajorticklabels):
            x,y = self.polar_transform(self.xymax,np.arccos(ccmajorticklabellocs[i]))
            self.ax.text(x,y,ccticklabel,color="black")
        for cctick in ccminorticks:
            r = [self.xymax-minorticklength,self.xymax            ]
            theta = [np.arccos(cctick),np.arccos(cctick)]
            x,y = self.polar_transform(r,theta)
            self.ax.plot(x,y,'-k')
 
        # re-draw xy ticks to be consistent with cc ticks (not elegant)
        for xx in xmajorticks:
            x = [xx,xx                         ]
            y = [0 ,majorticklength]
            self.ax.plot(x,y,'-k')
        for yy in ymajorticks:
            x = [0 ,majorticklength]
            y = [yy,yy                         ]
            self.ax.plot(x,y,'-k')


    def draw_point(self,ratio,cc,bias,bubblecolors=None,varlabels=None,labeloffset=0.025):
        
        if len(ratio.shape) == 2:
            nvars = len(ratio[:,0])
            ncases = len(ratio[0,:])
            for ivar in range(nvars):

                # transform coordinates
                r = ratio[ivar,:]
                theta = np.arccos(cc[ivar,:])
                size = np.abs(bias[ivar,:])/2.0
                x,y = self.polar_transform(r,theta)
                for icase in range(ncases):

                    # draw "bubbles" around points
                    if bubblecolors is None: 
                        bubblecolor='blue'
                    else:
                        bubblecolor=bubblecolors[icase]

                    circle = patches.Circle(
                            (x[icase],y[icase]),size[icase],
                            color=bubblecolor,
                            linewidth=0.0,
                            alpha=0.30,
                        )
                    self.ax.add_patch(circle)

                    # draw the actual points
                    pointcolor=bubblecolor
                    circle = patches.Circle(
                            (x[icase],y[icase]),0.01,
                            color=pointcolor,
                        )
                    self.ax.add_patch(circle)

                    # add labels to points if given
                    if varlabels is not None:
                        self.ax.text(
                                x[icase],y[icase]+labeloffset,
                                varlabels[ivar],
                                horizontalalignment='center',
                                verticalalignment='bottom',
                                color=bubblecolor,
                            )
        else:

            # transform coordinates
            r = ratio
            theta = np.arccos(cc)
            size = np.abs(bias)/2.0
            x,y = self.polar_transform(r,theta)
            for icase in range(len(x)):

                # draw "bubbles" around points
                if bubblecolors is None: 
                    bubblecolor='blue'
                else:
                    bubblecolor=bubblecolors[icase]

                circle = patches.Circle(
                        (x[icase],y[icase]),size[icase],
                        color=bubblecolor,
                        linewidth=0.0,
                        alpha=0.30,
                    )
                self.ax.add_patch(circle)

                # draw the actual points
                pointcolor=bubblecolor
                circle = patches.Circle(
                        (x[icase],y[icase]),0.01,
                        color=pointcolor,
                        #edgecolor="black",
                        #linewidth=0.5,
                    )
                self.ax.add_patch(circle)

                # add label to points if given
                if varlabels is not None:
                    self.ax.text(
                            x[icase],y[icase]+labeloffset,labels[icase],
                            horizontalalignment='center',
                            verticalalignment='bottom',
                            color=bubblecolor,
                        )

    def draw_rms(self,r,color='black',linestyle='dotted'):
        tmp = np.pi - np.arccos((r**2.0 + 1 - self.xymax**2.0)/2*r)
        theta1 = np.max([0.0,180.0/np.pi*tmp])
        theta2 = 180.0
        arc = patches.Arc(
                [1,0],2*r,2*r,
                theta1=theta1,
                theta2=theta2,
                color=color,
                linestyle=linestyle,
                linewidth=0.5,
            )
        self.ax.add_patch(arc)
