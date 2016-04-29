#!/usr/bin/env python3
'''
Copyright (C) 2015 Ryan Chiechi <r.c.chiechi@rug.nl>
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

import os,warnings,datetime
#TODO Replace csv with pandas
import csv
from shutil import copyfile
try:
    import numpy as np
    from scipy.interpolate import griddata
except ImportError as msg:
    pass #We catch numpy import errors in Parser.py
warnings.filterwarnings('ignore','.*comparison.*',FutureWarning)    

class Writer():
    '''The main Writer class for creating text files of parsed data.'''
    def __init__(self,parser):
        self.parser = parser
        self.opts = self.parser.opts
        if not os.path.exists(parser.opts.out_dir):
            print("Creating %s" % parser.opts.out_dir)
            os.mkdir(parser.opts.out_dir)

    def __getattr__(self, name):
        try:
                return getattr(self.parser, name) # 'inheret' the methods of self.parser
        except AttributeError as e:
                raise AttributeError("Writer object has no attribute '%s'" % name)

    def WriteParseInfo(self,extra=''):
        '''Write some summary information about the parameters
        used to parse the input data.'''
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_parseinfo.txt")
        with open(fn, 'a') as fh:
            fh.write("Parsed: %s\n" % str(datetime.datetime.today().ctime()) )
            t = str(vars(self.opts))
            t = t.replace(",","\n").replace("[","\n[")
            fh.write(t)
            fh.write(extra+"\n")    

    def WriteHistograms(self):
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Histograms.txt")
        with open(fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = []
            for x in self.X: headers += ["Log |J| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x]
            writer.writerow(headers)
            for i in range(0, len( self.XY[list(self.XY.keys())[0]]['hist']['bin'] ) ):
                row = []
                for x in self.X: row += ["%0.4f"%self.XY[x]['hist']['bin'][i], 
                                 "%s"%self.XY[x]['hist']['freq'][i],
                             "%0.4f"%self.XY[x]['hist']['fit'][i]]
                writer.writerow(row)
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_RHistograms.txt")
        with open(fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = []
            if self.opts.logr:
                for x in self.R['X']: headers += ["log |R| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x]
            else:
                for x in self.R['X']: headers += ["|R| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x]
            writer.writerow(headers)
            for i in range(0, len( self.R[list(self.R.keys())[0]]['hist']['bin'] ) ):
                row = []
                for x in self.R['X']: row += ["%0.4f"%self.R[x]['hist']['bin'][i],
                                 "%d"%self.R[x]['hist']['freq'][i],
                                 "%0.4f"%self.R[x]['hist']['fit'][i]]
                writer.writerow(row)

        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_LogdJdVHistograms.txt")
        with open(fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = []
            for x in self.X: headers += ["Log |dJ/dV| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x]
            writer.writerow(headers)
            for i in range(0, len( self.GHists[list(self.GHists.keys())[0]]['hist']['bin'] ) ):
                row = []
                for x in self.GHists: row += ["%0.4f"%self.GHists[x]['hist']['bin'][i], 
                                 "%s"%self.GHists[x]['hist']['freq'][i],
                             "%0.4f"%self.GHists[x]['hist']['fit'][i]]
                writer.writerow(row)


    def WriteVtrans(self):
        for key in ('pos', 'neg'):
            fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Vtrans_"+key+".txt")
            with open(fn, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, dialect='JV')
                writer.writerow(["Vtrans (eV)","Frequency",
                    "Gauss Fit (mean: %0.4f, Standard Deviation: %f)"%(self.FN[key]['mean'], self.FN[key]['std'])])
                data = {}
                for i in range(0, len(self.FN[key]['bin'])):
                    data[self.FN[key]['bin'][i]] = (self.FN[key]['freq'][i],self.FN[key]['fit'][i])
                for x in sorted(data.keys()):
                    writer.writerow(['%0.4f'%x,'%d'%data[x][0],'%0.2f'%data[x][1]])
                    
    def WriteFN(self):
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_FN.txt")
        with open(fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["1/V"] + ['Y_%d'%x for x in range(1,len( self.XY[list(self.XY.keys())[0]]['FN'] )+1)])
            for x in self.X:
                if x == 0.0:
                    continue
                writer.writerow(["%0.4f\t" % (1/x)] + list(self.XY[x]['FN']))
    
    def WriteGauss(self):
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Gauss.txt")
        with open(fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)","Log|J|","Standard Devaition","Standard Error of the Mean"])
            Y = []
            Yerr = []
            for x in self.X:
                writer.writerow(['%f'%x,'%f'%self.XY[x]['hist']['mean'],'%f'%self.XY[x]['hist']['std'],\
                        '%f'% (self.XY[x]['hist']['std']/np.sqrt(len(self.opts.in_files))) ])
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_RGauss.txt")
        with open(fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            if self.opts.logr:
                writer.writerow(["Potential (V)","log |R|","Standard Deviation"])
            else:
                writer.writerow(["Potential (V)","|R|","Standard Deviation"])
            Y = []
            Yerr = []
            for x in self.R['X']:
                writer.writerow(['%f'%x,'%f'%self.R[x]['hist']['mean'],'%f'%self.R[x]['hist']['std']])
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_logdJdVGauss.txt")
        with open(fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)","Log|dJ/dV|","Standard Devaition"])
            Y = []
            Yerr = []
            for x in self.GHists:
                writer.writerow(['%f'%x,'%f'%self.GHists[x]['hist']['mean'],'%f'%self.GHists[x]['hist']['std']])


    def WriteData(self, log=False):
        if log: key,label ='LogY','LogJ'
        else:   key, label ='Y','J'
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt")
        with open(fn,'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)"] + ['Y_%d'%x for x in range(1,len(self.XY[list(self.XY.keys())[0]][key] )+1)])
            for x in self.X:
                writer.writerow(["%0.4f"%x]+list(self.XY[x][key]))

    def WriteDJDV(self):
        label = 'DJDV'
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt")
        with open(fn,'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)"] + ['DJDV_%d'%x for x in range(1,len(self.DJDV[list(self.DJDV.keys())[0]])+1)])
            X = list(self.DJDV.keys())
            X.sort()
            for x in X:
                writer.writerow(["%0.4f"%x]+self.DJDV[x])

    def WriteFiltered(self):
        label = 'filtered'
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt")
        with open(fn,'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            for l in self.filtered:
                writer.writerow(l)

    def WriteRData(self):
        key,label = 'r', 'Rdata'
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt")
        with open(fn,'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)"] + ['Y_%d'%x for x in range(1,len(self.R[list(self.R.keys())[0]][key] )+1)])
            for x in self.R['X']:
                writer.writerow(["%0.4f"%x]+list(self.R[x][key]))

    def WriteGHistogram(self):
        '''Output for a contour plot of conductance'''
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_GHistogram.txt")
        with open(fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = ["Potential (V)", "Log dJ/dV", "Frequency"]
            writer.writerow(headers)

            for x in self.GHists:
                for i in range(0, len(self.GHists[x]['hist']['bin'])):
                    row = [x,"%0.4f"%self.GHists[x]['hist']['bin'][i],"%s"%self.GHists[x]['hist']['freq'][i]]
                    writer.writerow(row)
                writer.writerow([])

    def WriteGMatrix(self):
        '''Output for a matlab-style colormap maxtrix'''
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_GMatrix.txt")
        with open(fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            x,y,z = [],[],[]
            for i in range(0, len(self.GHists[list(self.GHists.keys())[0]]['hist']['bin'])):
                for v in self.GHists:
                    x.append(v)
                    y.append(self.GHists[v]['hist']['bin'][i])
                    z.append(self.GHists[v]['hist']['freq'][i])
            x,y,z = np.array(x),np.array(y),np.array(z)
            xmin,xmax = x.min(),x.max()
            ymin,ymax = self.opts.mlow, self.opts.mhi
            xi=np.linspace(xmin,xmax,200)
            yi=np.linspace(ymin,ymax,200)
            X,Y= np.meshgrid(xi,yi)
            Z = griddata((x, y), z, (X, Y),fill_value=0,method='cubic')
            
            #headers = ['MatrixData']
            #for x in X[0]:
            #   headers += ["%0.1f"%x]
            #writer.writerow(headers)
            
            for i in range(0,len(Z)):
                zi = []
                for z in Z[i]:
                    if z < 0:
                        zi.append(0)
                    else:
                        zi.append(z)
                #writer.writerow( ['%0.1f'%Y[i][0]]+list(Z[i]) )
                #writer.writerow( ['%0.1f'%Y[i][0]]+zi )
                writer.writerow(zi)
        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_GMatrix_Labels.txt")
        with open(fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = []
            for x in X[0]:
                headers += ["%0.1f"%x]
            writer.writerow(headers)
            headers = []
            for i in range(0,len(Y)):
                headers += ['%0.1f'%Y[i][0]]
            writer.writerow(headers)

    def WriteGeneric(self, dataset, bfn, labels=[]):
        '''Output for a generic set of data expecting an n-dimensional array'''

        if len(labels) and len(labels) != len(dataset):
            print("ERROR: Length of column labels does not match number of data columns for WriteGeneric!")
            return

        lencola = len(dataset[0])
        for d in dataset:
            if len(d) != lencola:
                print("ERROR: Length of columns differs for WriteGeneric!")
                return

        fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_"+bfn+".txt")
        with open(fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            if len(labels):
                writer.writerow(labels)

            for n in range(0, lencola):
                row = []
                for i in range(0,len(dataset)):
                    row.append(dataset[i][n])
                writer.writerow(row)


    def WriteGNUplot(self, gpinbn, tocopy=[]):
        absdir = os.path.dirname(os.path.abspath(__file__))
        tdir = os.path.join(absdir,'../templates/')
        gpintp = os.path.join(tdir,gpinbn+'.gpin')
        try:
            nsub = str(open(gpintp,'rt').read()).count('%s')
        except FileNotFoundError as msg:
            print("Could not read template file %s" % gpintp)
            return
        ssub = []
        for i in range(0,nsub):
            ssub.append(self.opts.outfile)
        txt = open(gpintp, 'rt').read() % tuple(ssub)
        fh = open(os.path.join(self.opts.out_dir,self.opts.outfile+'_'+gpinbn+'.gpin'), 'wt')
        fh.write(txt)
        fh.close()
        for fn in tocopy:
            copyfile(os.path.join(tdir,fn), \
                    os.path.join(self.opts.out_dir,fn))     
class Plotter():
    '''This is the main Plotter class for generating 
    plots using matplotlib.'''
    def __init__(self,parser):
        self.parser = parser
        self.opts = self.parser.opts

    def __getattr__(self, name):
        try:
                return getattr(self.parser, name)
        except AttributeError as e:
                raise AttributeError("Plotter object has no attribute '%s'" % name)

    def PlotData(self, key, ax, sym, **kw):
        xax = self.X
        if key == "FN":
            xax = 1/xax
            ax.set_title("Fowler Nordheim Plot of Initial Data")
            ax.set_xlabel(r'$V^{-1}$')
            ax.set_ylabel(r'$\mathregular{ln(\frac{J}{V^2})}$')
        if key == 'Y':
            if self.opts.compliance != np.inf: ax.set_ylim( (-1*self.opts.compliance, self.opts.compliance) )
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
                allY = np.append(allY,[self.XY[x][key][i] for x in self.X])
                ax.plot(xax,[self.XY[x][key][i] for x in self.X], sym, **kw)
            except IndexError:
                break
            except ValueError:
                break
        #if key == 'LogY':
        #   ax.set_ylim(allY.min(),allY.max())
        ax.axis([self.X.min(), self.X.max(), allY.min(),allY.max()])

    def PlotR(self, ax):
        if self.opts.logr: 
            ax.set_title("Semilog Plot of |R|")
            ax.set_ylabel(r'log|R|')
        else:
            ax.set_title("Plot of |R|")
            ax.set_ylabel(r'|R|')
        ax.set_xlabel("Potenial (V)")
        Y, Yerr = [],[]
        for x in self.R['X']:
            Y.append(self.R[x]["hist"]["mean"])
            Yerr.append(self.R[x]["hist"]["std"])
        ax.errorbar(self.R['X'], Y,  yerr=Yerr, marker='o', lw=0.0, color='k')

    def PlotRFit(self,ax):
        #print(self.R['X'])
        key = self.R['X'][-1]
        ax.set_title(r'Histogram and fit at $'+str(key)+'$ V')
        ax.set_xlabel(r'$log|R|$')
        ax.set_ylabel('Counts')
        ax.bar(self.R[key]['hist']['bin'], self.R[key]['hist']['freq'], width=0.05, color='r')
        ax.plot(self.R[key]['hist']['bin'], self.R[key]['hist']['fit'], lw=2.0, color='b', label='Fit')


    
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
        import matplotlib.pyplot as plt
        ax.set_title("Conductance Plot")
        ax.set_xlabel("Potential (V)")
        if self.opts.heatmapd == 0:
            ax.set_title("Heatmap of Initial Data")
            ax.set_ylabel(r'$log|J A cm^{-2}|$')
        elif self.opts.heatmapd == 1:
            ax.set_title("Derivative of Initial Data")
            ax.set_ylabel(r'$\mathregular{\frac{dJ}{dV}}$')
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
        ax.pcolormesh(X,Y,Z, cmap = plt.get_cmap('rainbow'))


    def PlotHist(self,ax):
        ax.set_title("Gaussian Fit and Raw Data")
        ax.set_xlabel('Potential (V)')
        ax.set_ylabel(r'Current Density $log_{10}|J(\mathrm{A cm}^{-2})|$')
        Y, Yerr = [],[]
        for x in self.X:
            Y.append(self.XY[x]["hist"]["mean"])
            Yerr.append(self.XY[x]["hist"]["std"])
        ax.errorbar(self.X, Y, yerr=Yerr, lw=3.0, color='k')

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


    def DoPlots(self, plt, *data):
        fig = plt.figure(figsize=(16,10))
        ax1 = fig.add_axes([0.06, 0.55, 0.4, 0.4])
        ax2 = fig.add_axes([0.56, 0.55, 0.4, 0.4])
        ax3 = fig.add_axes([0.06, 0.05, 0.4, 0.4])
        ax4 = fig.add_axes([0.56, 0.05, 0.4, 0.4])
        if 'J' in data:
            self.PlotFit(ax1)
            self.PlotData('LogY',ax2,':',lw=0.25, color='c')
            self.PlotHist(ax2)
        elif 'R' in data:
            self.PlotRFit(ax1)
            self.PlotR(ax2)
        self.PlotVtrans(ax4)
        self.PlotG(ax3)
        if self.opts.write:
            fig.savefig(os.path.join(self.opts.out_dir,self.opts.outfile+"_fig.png"), format="png")
        #self.PlotData('FN', ax3, 'x', ms=2)
        #self.PlotData('Y', ax1, '-')
        #self.PlotDJDV(ax1)


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


def WriteStats(out_dir, outfile, dataset, bfn, labels=[]):
    '''Output for a generic set of data expecting an n-dimensional array'''

    if len(labels) and len(labels) != len(dataset):
        print("ERROR: Length of column labels does not match number of data columns for WriteGeneric!")
        return

    lencola = len(dataset[0])
    for d in dataset:
        if len(d) != lencola:
            print("ERROR: Length of columns differs for WriteGeneric!")
            return
    fn = os.path.join(out_dir,outfile+"_"+bfn+".txt")
    with open(fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        if len(labels):
            writer.writerow(labels)

        for n in range(0, lencola):
            row = []
            for i in range(0,len(dataset)):
                row.append(dataset[i][n])
            writer.writerow(row)

