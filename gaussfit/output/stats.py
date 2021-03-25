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
import csv
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
    # from scipy.interpolate import griddata
except ImportError as msg:
    pass #We catch numpy import errors in Parser.py
warnings.filterwarnings('ignore','.*comparison.*',FutureWarning)


class StatPlotter:
    '''This is the main Plotter class for Stats.py.'''
    def __init__(self, statparser):
        self.opts = statparser.opts
        self.dataset = statparser.dataset
        self.cutoff = 0.01

    def PlotGmeanData(self, key, ax, sym, **kw):

        xmin = np.asarray(self.dataset[key][0]).min()
        xmax = np.asarray(self.dataset[key][0]).max()
        ymin = np.asarray(self.dataset[key][1]).min()
        ymax = np.asarray(self.dataset[key][1]).max()
        if ymax < self.cutoff:
            ymax = self.cutoff
        if ymin > self.cutoff:
            ymin = self.cutoff
        ymax = ymax + ymax*0.1
        ymin = ymin - ymin*0.1
        ax.set_yscale('log')
        ax.set_title("P-values for %s (Gmean)" % key)
        ax.set_xlabel("Potenial (V)")
        ax.set_ylabel('p-value')

        ax.plot(self.dataset[key][0],self.dataset[key][1], sym, **kw)
        ax.axis([xmin,xmax,ymin,ymax])

        self.PlotCutoff(xmin, xmax, ax)

    def PlotJData(self, key, ax, sym, **kw):

        xmin = np.asarray(self.dataset[key][0]).min()
        xmax = np.asarray(self.dataset[key][0]).max()
        ymin = np.asarray(self.dataset[key][3]).min()
        ymax = np.asarray(self.dataset[key][3]).max()
        if ymax < self.cutoff:
            ymax = self.cutoff
        if ymin > self.cutoff:
            ymin = self.cutoff
        ymax = ymax + ymax*0.1
        ymin = ymin - ymin*0.1
        ax.set_yscale('log')
        ax.set_title("P-values for %s (mean using N)" % key)
        ax.set_xlabel("Potenial (V)")
        ax.set_ylabel('p-value')

        ax.plot(self.dataset[key][0],self.dataset[key][3], sym, **kw)
        ax.axis([xmin,xmax,ymin,ymax])

        self.PlotCutoff(xmin, xmax, ax)

    def PlotCutoff(self, xmin, xmax,  ax):
        X = np.linspace(xmin,xmax)
        ax.plot(X, [ self.cutoff for x in X],  '-', lw='0.5', color='c')

    def DoPlots(self, plt):
        fig = plt.figure(figsize=(10,5))
        #ax1 = fig.add_axes([0.08, 0.1, 0.4, 0.8])
        #ax2 = fig.add_axes([0.56, 0.1, 0.4, 0.8])
        ax1 = fig.add_axes([0.06, 0.55, 0.4, 0.35])
        ax2 = fig.add_axes([0.56, 0.55, 0.4, 0.35])
        ax3 = fig.add_axes([0.06, 0.05, 0.4, 0.35])
        ax4 = fig.add_axes([0.56, 0.05, 0.4, 0.35])
        self.PlotGmeanData('J',ax1,'.',lw=1.25, color='b')
        self.PlotGmeanData('R',ax2,'.',lw=1.25, color='b')
        self.PlotJData('J',ax3,'.',lw=1.25, color='b')
        self.PlotJData('R',ax4,'.',lw=1.25, color='b')
        if self.opts.write:
            fig.savefig(self.opts.outfile+"_statfig.png", format="png")






def WriteStats(out_dir, outfile, dataset, bfn, labels=None):
    '''Output for a generic set of data expecting an n-dimensional array'''

    labels = labels or []

    if len(labels) and len(labels) != len(dataset):
        logger.error("Length of column labels does not match number of data columns for WriteGeneric!")
        return

    lencola = len(dataset[0])
    for d in dataset:
        if len(d) != lencola:
            logger.error("Length of columns differs for WriteGeneric!")
            return
    _fn = os.path.join(out_dir,outfile+"_"+bfn+".txt")
    with open(_fn , 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        if len(labels):
            writer.writerow(labels)

        for n in range(0, lencola):
            row = []
            for i in range(0,len(dataset)):
                row.append(dataset[i][n])
            writer.writerow(row)