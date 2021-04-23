'''
Copyright (C) 2018 Ryan Chiechi <r.c.chiechi@rug.nl>
Description:

    This absolute mess outputs parsed data into text files and plots. If
    you make changes, do so with care!

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import os
import warnings
# import datetime
#TODO Replace csv with pandas
# import csv
# from shutil import copyfile
import logging
# from scipy.special import stdtrit #pylint: disable=E0611
# import matplotlib.pyplot as plt
from gaussfit.colors import GREEN, TEAL, YELLOW, WHITE


logger = logging.getLogger('output')
loghandler = logging.StreamHandler()
loghandler.setFormatter(logging.Formatter(\
                fmt=GREEN+os.path.basename('%(name)s'+TEAL)+' %(levelname)s '+YELLOW+'%(message)s'+WHITE))
logger.addHandler(loghandler)

try:
    import numpy as np
    from scipy.interpolate import griddata
    from scipy.special import stdtrit #pylint: disable=E0611
    from matplotlib import container
except ImportError as msg:
    pass #We catch numpy import errors in Parser.py
warnings.filterwarnings('ignore','.*comparison.*',FutureWarning)




class Plotter():
    '''
    This is the main Plotter class for generating
    plots using matplotlib.
    '''

    def __init__(self,parser,plt):
        self.parser = parser
        self.opts = self.parser.opts
        self.plt = plt

    def __getattr__(self, name):
        try:
            return getattr(self.parser, name)
        except AttributeError as msg:
            raise AttributeError("Plotter object has no attribute '%s'" % name) from msg

    def PlotData(self, key, ax, sym, **kw):
        xax = np.array(list(self.XY))
        if key == "FN":
            xax = 1/xax
            ax.set_title("Fowler Nordheim Plot of Initial Data")
            ax.set_xlabel(r'$V^{-1}$')
            ax.set_ylabel(r'$\mathregular{ln(\frac{J}{V^2})}$')
        if key == 'Y':
            if self.opts.compliance != np.inf:
                ax.set_ylim( (-1*self.opts.compliance, self.opts.compliance) )
            ax.set_title("Initial Data")
            ax.set_xlabel("Potenial (V)")
            ax.set_ylabel(r'Current Density ($A cm^{-2}$)')
        if key == 'LogY':
            ax.set_title("Semilog Plot of Initial Data")
            ax.set_xlabel("Potenial (V)")
            ax.set_ylabel(r'Current Density $\mathregular{log_{10}|J(\mathrm{A cm^{-2}})|}$')
        i = -1
        allY = np.array([])
        while True:
            i += 1
            try:
                allY = np.append(allY,[self.XY[x][key][i] for x in self.XY])
                ax.plot(xax,[self.XY[x][key][i] for x in self.XY], sym, **kw)
            except IndexError:
                break
            except ValueError:
                break
        #if key == 'LogY':
        #   ax.set_ylim(allY.min(),allY.max())
        ax.axis([xax.min(), xax.max(), allY.min(),allY.max()])

    def PlotSegmentedGauss(self, ax):
        '''Plot segmented J/V data'''
        ax.set_title("Semilog Plots of |J| by Trace")
        ax.set_ylabel(r'Current Density $log_{10}|J(\mathrm{A cm}^{-2})|$')
        ax.set_xlabel(r'Potential (V)')
        X,Y={},{}
        Yerr={}
        for segment in self.segments:
            for trace in self.segments[segment]:
                #TODO: Fix this hack
                if trace == 'combined':
                    continue
                if trace not in X:
                    X[trace] = []
                    Y[trace] = []
                    Yerr[trace] = []
                # _x = list(self.segments[segment][trace].keys())
                # _x.sort()
                for x in self.segments[segment][trace]:
                    X[trace].append(x)
                    _hist = self.segments[segment][trace][x]
                    Y[trace].append(_hist['mean'])
                    _sem = float(_hist['std'])/np.sqrt(len(self.opts.in_files)-1 or 1)
                    _t_val = _sem * stdtrit( len(self.opts.in_files)-1 or 1, 1 - self.opts.alpha )
                    Yerr[trace].append(_t_val)
        for trace in Y:
            ax.errorbar(X[trace],
                        Y[trace],
                        yerr=Yerr[trace],
                        marker='o',
                        lw=0,
                        elinewidth=0.25,
                        capsize=0.5,
                        label='Trace %s' % (trace+1))
        handles, labels = ax.get_legend_handles_labels()
        handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
        ax.legend(handles, labels)
    def PlotR(self, ax):
        if self.opts.logr:
            ax.set_title("Semilog Plot of |R|")
            ax.set_ylabel(r'log|R|')
        else:
            ax.set_title("Plot of |R|")
            ax.set_ylabel(r'|R|')
        ax.set_xlabel("Potenial (V)")
        Y, Yerr = [],[]
        for x in self.XY:
            Y.append(self.XY[x]['R']["hist"]["mean"])
            _sem = float(self.XY[x]['R']["hist"]["std"])/np.sqrt(len(self.opts.in_files)-1 or 1)
            _t_val = _sem * stdtrit( len(self.opts.in_files)-1 or 1, 1 - self.opts.alpha )
            Yerr.append(_t_val)
        ax.errorbar(list(self.XY), Y,  yerr=Yerr, marker='o', lw=0.0, color='k')

    def PlotRFit(self,ax):
        key = list(self.XY)[-1]
        ax.set_title(r'Histogram and fit at $'+str(key)+'$ V')
        ax.set_xlabel(r'$log|R|$')
        ax.set_ylabel('Counts')
        ax.bar(self.XY[key]['R']['hist']['bin'], self.XY[key]['R']['hist']['freq'], width=0.05, color='r')
        ax.plot(self.XY[key]['R']['hist']['bin'], self.XY[key]['R']['hist']['fit'], lw=2.0, color='b', label='Fit')



    def PlotDJDV(self,ax):
        xax = list(self.DJDV.keys())
        xax.sort()
        ax.set_xlabel("Potential (V)")
        ax.set_title("Derivative of Initial Data")
        ax.set_ylabel(r'$\mathregular{\frac{dJ}{dV}}$')
        #ax.set_ylabel(r'Normalized $\mathregular{\frac{dJ}{dV}}$')
        #ax.axis([np.array(xax).min(), np.array(xax).max(), pow(10,self.opts.mlow), pow(10,self.opts.mhi)])
        i = -1
        while True:
            i += 1
            try:
                if i not in self.ohmic:
                    ax.plot(xax, [self.DJDV[x][i] for x in xax], "-", lw=2)
                else:
                    ax.plot(xax, [self.DJDV[x][i] for x in xax], "-", lw=0.5, color='grey')
            except IndexError:
                break

    def PlotG(self,ax):
        # import matplotlib.pyplot as plt
        ax.set_title("Conductance Plot")
        ax.set_xlabel("Potential (V)")
        if self.opts.heatmapd == 0:
            ax.set_title("Heatmap of Initial Data")
            ax.set_ylabel(r'$log|J A cm^{-2}|$')
        elif self.opts.heatmapd == 1:
            ax.set_title("Derivative of Initial Data")
            ax.set_ylabel(r'$\mathregular{\log|\frac{dJ}{dV}|}$')
        elif self.opts.heatmapd == 2:
            ax.set_title("Second Derivative of Initial Data")
            ax.set_ylabel(r'$\mathregular{\frac{d^2J}{dV^2}}$')
        #ax.set_ylabel(r'$log|\mathregular{\frac{dJ}{dV}}|$')
        x,y,z =[],[],[]
        for v in self.GHists:
            for i in range(0, len(self.GHists[v]['hist']['bin'])):
                if i in self.ohmic and self.opts.skipohmic:
                    continue
                x.append(v)
                y.append(self.GHists[v]['hist']['bin'][i])
                z.append(self.GHists[v]['hist']['freq'][i])
        x,y,z = np.array(x),np.array(y),np.array(z)
        xmin,xmax = x.min(),x.max()
        #ymin,ymax = y.min(),y.max()
        ymin,ymax = self.opts.mlow, self.opts.mhi
        x = np.r_[x,xmin,xmax]
        y = np.r_[y,ymin,ymax]
        z = np.r_[z,z[0],z[-1]]
        xi = np.linspace(xmin, xmax, 200)
        yi = np.linspace(ymin, ymax, 200)
        X,Y= np.meshgrid(xi,yi)
        Z = griddata((x, y), z, (X, Y),method='nearest')
        ax.axis([xmin, xmax, ymin, ymax])
        ax.pcolormesh(X,Y,Z, cmap = self.plt.get_cmap('rainbow'))


    def PlotNDC(self,ax):
        # import matplotlib.pyplot as plt
        ax.set_title("NDC Plot")
        ax.set_xlabel("Potential (V)")
        ax.set_title("Heatmap of NDC")
        ax.set_ylabel('Normalized Differential Condutance')
        x,y,z =[],[],[]
        for v in self.NDCHists:
            for i in range(0, len(self.NDCHists[v]['hist']['bin'])):
                if i in self.ohmic and self.opts.skipohmic:
                    continue
                x.append(v)
                y.append(self.NDCHists[v]['hist']['bin'][i])
                z.append(self.NDCHists[v]['hist']['freq'][i])
        x,y,z = np.array(x),np.array(y),np.array(z)
        xmin,xmax = x.min(),x.max()
        ymin,ymax = self.opts.ndc_mlow, self.opts.ndc_mhi
        x = np.r_[x,xmin,xmax]
        y = np.r_[y,ymin,ymax]
        z = np.r_[z,z[0],z[-1]]
        xi = np.linspace(xmin, xmax, 200)
        yi = np.linspace(ymin, ymax, 200)
        X,Y= np.meshgrid(xi,yi)
        Z = griddata((x, y), z, (X, Y),method='nearest')
        ax.axis([xmin, xmax, ymin, ymax])
        ax.pcolormesh(X,Y,Z, cmap = self.plt.get_cmap('rainbow'))


    def PlotHist(self,ax):
        ax.set_title("Gaussian Fit and Raw Data")
        ax.set_xlabel('Potential (V)')
        ax.set_ylabel(r'Current Density $log_{10}|J(\mathrm{A cm}^{-2})|$')
        Y, Yerr = [],[]
        for x in self.XY:
            Y.append(self.XY[x]["hist"]["mean"])
            Yerr.append(self.XY[x]["hist"]["std"])
        ax.errorbar(list(self.XY.keys()), Y, yerr=Yerr, lw=3.0, color='k')

    def PlotVtrans(self,ax):
        ax.set_title(r'Histogram and fit of $V_{trans}$')
        ax.set_xlabel(r'$V_{trans}$')
        ax.set_ylabel('Counts')
        for key in ('pos','neg'):
            ax.bar(self.FN[key]['bin'], self.FN[key]['freq'], width=0.01, color='g')
            ax.plot(self.FN[key]['bin'], self.FN[key]['fit'], lw=2.0, color='b', label='Fit')

    def PlotFit(self,ax):
        key = list(self.XY.keys())[-1]
        ax.set_title(r'Histogram and fit at $'+str(key)+'$ V')
        ax.set_xlabel(r'$log|J|$')
        ax.set_ylabel('Counts')
        ax.bar(self.XY[key]['hist']['bin'], self.XY[key]['hist']['freq'], width=0.1, color='r')
        ax.plot(self.XY[key]['hist']['bin'], self.XY[key]['hist']['fit'], lw=2.0, color='b', label='Fit')


    def DoPlots(self):
        fig = self.plt.figure(figsize=(16,10))
        ax1 = fig.add_axes([0.06, 0.55, 0.4, 0.4])
        ax2 = fig.add_axes([0.56, 0.55, 0.4, 0.4])
        ax3 = fig.add_axes([0.06, 0.05, 0.4, 0.4])
        ax4 = fig.add_axes([0.56, 0.05, 0.4, 0.4])
        if self.opts.plots == 'J':
            self.PlotFit(ax1)
            self.PlotData('LogY',ax2,':',lw=0.25, color='c')
            self.PlotHist(ax2)
        elif self.opts.plots == 'R':
            self.PlotRFit(ax1)
            self.PlotR(ax2)
        #self.PlotVtrans(ax4)
        self.PlotSegmentedGauss(ax4)
        if self.opts.histplots == 'NDC':
            self.PlotNDC(ax3)
        elif self.opts.histplots == 'G':
            self.PlotG(ax3)
        if self.opts.write:
            fig.savefig(os.path.join(self.opts.out_dir,self.opts.outfile+"_fig.png"), format="png")
        #self.PlotData('FN', ax3, 'x', ms=2)
        #self.PlotData('Y', ax1, '-')
        #self.PlotDJDV(ax1)
