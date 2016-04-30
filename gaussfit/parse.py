#!/usr/bin/env python3
'''
Copyright (C) 2016 Ryan Chiechi <r.c.chiechi@rug.nl>
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

import sys,os,logging,warnings,csv,threading
from gaussfit.colors import *
from gaussfit.logger import DelayedHandler
#import concurrent.futures 

try:
    import pandas as pd
    from scipy.optimize import curve_fit,OptimizeWarning
    import scipy.interpolate 
    from scipy.stats import gmean,norm
    import scipy.misc 
    import numpy as np
    # SciPy throws a useless warning for noisy J/V traces
    warnings.filterwarnings('ignore','.*Covariance of the parameters.*',OptimizeWarning)

except ImportError as msg:
    print("\n\t\t%s> > > Error importing numpy/pandas/scipy! %s%s%s < < <%s" % (RED,RS,str(msg),RED,RS))
    sys.exit()

warnings.filterwarnings('ignore','.*divide by zero.*',RuntimeWarning)
warnings.filterwarnings('ignore','.*',UserWarning)
#warnings.filterwarnings('ignore','.*invalid value encountered in log.*',RuntimeWarning)
#warnings.filterwarnings('ignore','.*invalid value encountered in true_divide.*',RuntimeWarning)
warnings.filterwarnings('ignore','.*invalid value encountered.*',RuntimeWarning)
warnings.filterwarnings('ignore','.*Mean of empty slice.*', RuntimeWarning)
#warnings.filterwarnings('error','.*Degrees of freedom <= 0 for slice.*', RuntimeWarning)
#warnings.filterwarnings('ignore','.*impossible result.*',UserWarning)




class Parse():
    '''
    This is the main parsing class that takes input data
    in the form of text files, parses them into dictionary
    and list attributes and provides methods for performing
    operations: Gaussian fits, F-N calculations, Rectification,
    Vtrans, and dJ/DV.
    '''
    def __init__(self,opts,handler=None,lock=None):
        self.error = False
        self.opts = opts
        self.df = pd.DataFrame()
        self.XY = {}
        self.X = np.array([])
        self.FN = {}
        self.compliance_traces = []
        self.ohmic = []
        self.DJDV = {}
        self.GHists = {}
        self.filtered = []
        self.R = {}
        self.traces = {}
        # Pass a lock when calling inside a thread
        if lock:
            self.lock = lock
        else:
            self.lock = threading.Lock()
        # Pass your own log hanlder, e.g., when calling from a GUI
        if not handler:
            self.loghandler = DelayedHandler()
            self.loghandler.setFormatter(logging.Formatter(\
                fmt=GREEN+os.path.basename(sys.argv[0]+TEAL)+' %(levelname)s '+YELLOW+'%(message)s'+WHITE))
        else:
            self.loghandler = handler
        self.logger = logging.getLogger('parser')
        self.logger.addHandler(self.loghandler)
        self.logger.setLevel(getattr(logging,self.opts.loglevel.upper()))

        if self.opts.nomin:
            self.logger.debug("Using interpolation on FN")

    @classmethod
    def doOutput(cls,writer):
        ''' Run through all the methods for writing output files.'''
        writer.logger.info("Writing files...")
        writer.WriteParseInfo()
        writer.WriteVtrans()
        writer.WriteGNUplot('Vtransplot')
        writer.WriteFN()
        writer.WriteGauss()
        writer.WriteGNUplot('JVplot')
        writer.WriteData()
        writer.WriteDJDV()
        writer.WriteFiltered()
        writer.WriteData(True)
        writer.WriteRData()
        writer.WriteGNUplot('Rplot')
        try:
            writer.WriteHistograms()
            writer.WriteGNUplot('JVhistplot')
        except IndexError:
            print("Error outputting histrograms")
        try:
            writer.WriteGHistogram()
        except IndexError:
            print("Error outputting Ghistrograms")
        try:
            writer.WriteGMatrix()
            writer.WriteGNUplot('Gplot', ['parula.pal']) 
        except IndexError:
            print("Error outputting GMatrix")
        writer.logger.info("Done!")

    def ReadFiles(self, fns):
        '''Walk through input files and parse
        them into attributes '''
        if type(fns) == type(str()):
            fns = [fns]
        if self.opts.Ycol > 0:
            self.logger.info("Parsing two columns of data.")
            self.logger.debug('Parsing %s' % ', '.join(fns))
            frames = {}
            for f in fns:
                try:
                    frames[f]=pd.read_csv(f,sep=self.opts.delim,usecols=(self.opts.Xcol,self.opts.Ycol),names=('V','J'),header=0)
                except OSError as msg:
                    self.logger.warn("Skipping %s because %s" % (f,str(msg)))
            self.df = pd.concat(frames)
        else:
            #TODO This needs to be multiindex
            self.logger.info("Parsing all columns of data.")
            self.df = pd.concat((pd.read_csv(f,sep=self.opts.delim,
                header=None,skiprows=1) for f in fns),ignore_index=True)
            X,Y = [],[]
            for row in self.df.iterrows():
                for y in row[1][1:]:
                    X.append(row[1][0])
                    Y.append(y)
            self.df = pd.DataFrame({'V':X,'J':Y})
        self.__parse()

    def ReadPandas(self,df):
        '''Take a pandas.DataFrame as input instead of files.'''
        self.logger.debug("Using Pandas as input")
        self.df = df
        self.__parse()

    def __parse(self):
        '''Read a pandas.DataFrame and compute Fowler-Nordheim
        values, log10 the J or I values and create a dictionary
        indexed by unique voltages.'''
        if (self.df.V.dtype,self.df.J.dtype) != ('float64','float64'):
            self.logger.error("Parsed data does not appear to contain numerical data!")
            self.error = True
            return
        try:
            self.df['FN'] = np.log(abs(self.df.J)/self.df.V**2)
        except ZeroDivisionError:
            self.logger.warn("Error computing FN (check your input data).")
            self.df['FN'] = np.array([x*0 for x in range(0, len(self.df['V']))])
        self.df.J.replace(0.0,value=10e-16,inplace=True)
        self.df['logJ'] = np.log10(abs(self.df.J)) # Cannot log10 zero
        self.logger.info('%s values of log|J| above compliance (%s)' % 
                (len(self.df['logJ'][self.df['logJ']>self.opts.compliance]),self.opts.compliance))
        
        #The default log handler only emits when you call flush() after setDelay() called
        self.loghandler.setDelay()
        
        for x, group in self.df.groupby('V'):
            self.XY[x] = { "Y":group['J'], 
                   "LogY":group['logJ'], 
                   "hist":self.dohistogram(group['logJ'],"J"), 
                   "FN": group['FN']}
        
        self.X = sorted(self.XY.keys())
        self.logger.debug("X = %s" % str(self.X) ) 
        self.X = np.array(self.X)
        
        self.logger.info("Done parsing input data")
        self.logger.info("* * * * * * Finding traces   * * * * * * * *")
        self.loghandler.flush()
        self.findTraces()
        self.logger.info("* * * * * * Computing dY/dX  * * * * * * * *")
        self.loghandler.flush()
        #NOTE Threads do not seem to help with performance
        t1 = threading.Thread(target = self.dodjdv) # This must come first for self.ohmic to be populated!
        t1.start()
        if self.opts.skipohmic:
                t1.join()
        self.logger.info("* * * * * * Computing Vtrans * * * * * * * *")
        self.loghandler.flush()
        t2 = threading.Thread(target=self.findmin)
        t2.start()
        self.logger.info("* * * * * * Computing |R|  * * * * * * * * *")
        self.loghandler.flush()
        t3 = threading.Thread(target=self.dorect)
        t3.start()
        t1.join()
        t2.join()
        t3.join()
        self.logger.info("* * * * * * * * * * * * * * * * * * * * * * ")
        self.loghandler.flush()
        self.PrintFN()
        self.loghandler.unsetDelay()

    def findTraces(self):
        '''Try to find individual J/V traces. A trace is defined
        by a complete forward and reverse trace unless the input
        dataset comprises only forward or only reverse traces.
        findTraces will average the two sweeps, which may cause
        problems with very hysteretic data.'''
        def __checktraces(traces):
            if not traces:
                self.logger.error("No traces!?") # This should never happen!
                return False
            v = self.df.loc[self.df.index.levels[0][0]]['V'].values
            st,ed = [v[0]],[v[-1]]
            for r in self.df.index.levels[0][1:]:
                s,e = self.df.loc[r]['V'].values[0],self.df.loc[r]['V'].values[-1]
                if s not in st or e not in ed:
                    self.logger.warn('file "%s" starts and ends with weird voltages (%s -> %s)' % (r,s,e))
                st.append(s)
                ed.append(e)
            if self.opts.tracebyfile:
                return True
            self.logger.debug("Checking V starting from slice %s:%s" % (traces[0][0],traces[0][1]) )
            lt = len(self.df.V[ traces[0][0]:traces[0][1] ])
            for trace in traces:
                if lt != len(self.df.V[trace[0]:trace[1]]):
                    #TODO Figure out how to mask/delete these files from parsing
                    if trace[0][0] != trace[1][0]:
                        self.logger.warn('Unequal voltage steps somewhere bewteen "%s" (and) "%s"' % (trace[0][0],trace[1][0]) )
                    else:
                        self.logger.warn('Unequal voltage steps somewhere in "%s"' % trace[0][0] )
                    self.loghandler.flush()
                    return False
                self.logger.debug("Trace: %s -> %s" % (self.df.V[trace[0]],self.df.V[trace[-1]]) )
            self.logger.info("Traces look good.")
            return True
        self.loghandler.flush() 
        traces = []
        ntraces = 0

        if self.opts.tracebyfile:
            self.logger.info("Assuming each file contains one (foward/backward) trace.")
            for r in self.df.index.levels[0]:
                traces.append( ((r,self.df.loc[r].index[0]),(r,self.df.loc[r].index[-1])) )
            ntraces = len(traces)

        if not ntraces and self.df.V.value_counts().index[0] == 0.0:
            #NOTE t is a tuple with both indices 0 = filename, 1 = index
            try:
                ntraces = int(self.df.V.value_counts()[0]/3) # Three zeros in every trace!
                self.logger.info("This looks like an EGaIn dataset.")
                for t in zip(*(iter(self.df[self.df.V == 0.00].V.index),) * 3):
                    traces.append( (t[0],t[2]) )
            except ValueError:
                self.logger.warn("Did not find three-zero (EGaIn) traces!")
        if not ntraces:
            self.logger.warn("This does not look like an EGaIn dataset.")
            try:
                ntraces = int(self.df.V.value_counts()[self.df.V[0]]/2) # Two end-pointss in every trace!
                for t in zip(*(iter(self.df[self.df.V == self.df.V[0]].V.index),) * 2):
                    traces.append( (t[0],t[1]) )
            except ValueError:
                self.logger.warn("Did not find three-zero (EGaIn) traces!")

        if not ntraces or not __checktraces(traces):
            self.logger.warn("Recomputing traces based on repeat values")
            traces = []
            Vinit = self.df.V[0]
            trace = []
            for row in self.df[0:].iterrows():
                if row[1].V == Vinit:
                    if not len(trace) or len(trace) == 1:
                        trace.append(row[0])
                    elif len(trace) == 2:
                        traces.append(trace)
                        #print(self.df[trace[0]:trace[1]+1])
                        trace = [row[0]]
            ntraces = len(traces)
        if not __checktraces(traces):
            self.logger.error("Problem with traces: FN and derivative probably will not work correctly!")
            self.loghandler.flush()
            #for r in self.df.index.levels[0]:
            #    self.logger.info('(%s) start: %s end: %s' % (r,
            #            self.df.loc[r]['V'].values[0],self.df.loc[r]['V'].values[-1]))
            #self.loghandler.flush()
        self.logger.info("Found %s traces (%s)." % (ntraces,len(traces)) )
        self.traces = traces
        self.loghandler.flush() 
        idx = []
        frames = {}
        self.logger.info("Compressing forward/reverse sweeps to single traces.")
        self.loghandler.flush()
        for col in range(0,len(self.traces)):
            fbtrace = self.df[self.traces[col][0]:self.traces[col][1]].sort_values('V')
            avg = {'V':[],'J':[],'FN':[]}
            for x,group in fbtrace.groupby('V'):
                avg['V'].append(x)
                avg['J'].append(self.signedgmean(group['J']))
                #avg['J'].append(np.mean(group['J']))
                if not self.opts.nomin:
                    fn = np.mean(group['FN'][group['FN'] <= self.opts.compliance])
                    #NOTE signedgmean was very slow at one point so I disabled it
                    #fn = self.signedgmean(group['FN'][group['FN'] <= self.opts.compliance])
                else:
                    fn = np.mean(group['FN'])
                    #fn = self.signedgmean(group['FN'])
                avg['FN'].append(fn)
            frames[col] = pd.DataFrame(avg)
        try:
            self.avg = pd.concat(frames)
        except ValueError:
            self.error = True
            self.avg = pd.DataFrame()
            self.logger.error('Unable to parse traces.')
    def dodjdv(self):
        '''
        Fit a spline function to X/Y data, 
        compute dY/dX and normalize.
        '''
        linx = np.linspace(self.df.V.min(), self.df.V.max(), 200)
        if self.opts.vcutoff > 0:
            vfilterneg,vfilterpos = np.linspace(-1*self.opts.vcutoff,0,200), np.linspace(0,self.opts.vcutoff.max(),200)
        else:
            vfilterneg,vfilterpos = np.linspace(self.df.V.min(),0,200), np.linspace(0,self.df.V.max(),200)
        if self.opts.vcutoff > 0:
            vfilterneg,vfilterpos = linx[-1*self.opts.vcutoff < linx < 0],linux[0 < linx < self.opts.vcutoff]
        else:
            vfilterneg,vfilterpos = linx[linx < 0], linx[linx > 0]

        spls = {}
        splhists = {}
        filtered = [('Potential', 'Fit', 'Y')]
        for x in linx: 
            spls[x] = []
            splhists[x] = {'spl':[],'hist':{}}
        
        for trace in self.avg.index.levels[0]:
            try:
                spl = scipy.interpolate.UnivariateSpline(self.avg.loc[trace]['V'],self.avg.loc[trace]['J'], k=5, s=self.opts.smooth )
                dd =  scipy.interpolate.UnivariateSpline(self.avg.loc[trace]['V'],self.avg.loc[trace]['J'], k=5, s=None).derivative(2)
            except Exception as msg:
                self.logger.error('Error in derivative calulation: %s' % str(msg))
                continue

            spldd = dd(vfilterpos) #Compute d2J/dV2
            spldd += -1*dd(vfilterneg) #Compute d2J/dV2
            if len(spldd[spldd<0]):
                    # record in the index where dY/dX is < 0 within vcutoff range
                    self.ohmic.append(trace)  
                    if self.opts.skipohmic:
                        continue
            else:
                for row in self.avg.loc[trace].iterrows():
                    # filtered is a list containing only "clean" traces         
                    filtered.append( (row[1].V, spl(row[1].V), row[1].J) )
            err = None
            for x in sorted(spls.keys()):
                try:
                    d = spl.derivatives(x)
                except ValueError as msg:
                    err = str(msg)
                    self.logger.warn('Error computing derivative: %s' % str(msg))
                    continue
                if np.isnan(d[self.opts.heatmapd]):
                    self.logger.warn("Got NaN computing dJ/dV")
                    continue
                spls[x].append(d[self.opts.heatmapd])
                splhists[x]['spl'].append(np.log10(abs(d[self.opts.heatmapd])))
            if err:
                self.logger.error("Error while computing derivative: %s" % str(err))

        self.logger.info("Non-tunneling traces: %s (out of %0d)" % 
                    ( len(self.ohmic), len(self.traces) ) )
        self.loghandler.flush()
        for x in splhists:
            splhists[x]['hist'] = self.dohistogram(np.array(splhists[x]['spl']), label='DJDV')
        self.logger.info("dJdV complete.")
        self.loghandler.flush()
        self.DJDV, self.GHists, self.filtered = spls, splhists, filtered
        sys.exit()

    def dorect(self):
        ''' 
        Divide each value of Y at +V by Y at -V
        and build a histogram of rectification, R
        '''
        R = {}
        r = {}
        clipped = 0
        for trace in self.traces:
            rows = {}
            for row in self.df[trace[0]:trace[1]].iterrows():
                if row[1].V in rows:
                    #TODO Don't just average out hysterysis, it is important for R
                    rows[row[1].V] = np.mean([rows[row[1].V],row[1].J])
                else:
                    rows[row[1].V] = row[1].J
                if row[1].V not in r:
                    r[row[1].V] = []
            for x in rows:
                #TODO AFM data break this!
                if x not in rows or -1*x not in rows:
                    self.logger.warn("Rectification data missing voltages.")
                    continue
                if x == 0.0:
                    r[x].append(1.)
                    continue
                if self.opts.logr:
                    r[x].append(np.log10(abs(rows[x]/rows[-1*x])))
                else:
                    r[x].append(abs(rows[x]/rows[-1*x]))
                if r[x][-1] > self.opts.maxr:
                    clipped += 1
        for x in r:
            if x < 0: continue
            R[x]={'r':np.array(r[x]),'hist':self.dohistogram(np.array(r[x]),"R")}
            R[-1*x] = R[x]
        R['X'] = np.array(sorted(R.keys()))
        if clipped:
            if self.opts.logr: rstr = 'log|R|'
            else: rstr = '|R|'
            self.logger.info("%s values of %s exceed maxR (%s)" % (clipped, rstr, self.opts.maxr))
        self.logger.info("R complete.")
        self.loghandler.flush()
        self.R = R
    
    def getminroot(self, spl):
        '''
        Return the root of the first derivative of a spline function
        that results in the lowest value of ln(Y/X^2) vs. 1/X
        '''
        splvals = {}
        for r in spl.derivative().roots():
            splvals[float(spl(r))] = r
        if not splvals:
            return np.NAN
            #return False
        # Because UnivariateSpline crashes when we use 
        # 1/X, the plots are flipped so we take the max
        return splvals[np.nanmax(list(splvals.keys()))]

    def findmin(self):
        '''
        Find the troughs of ln(Y/X^2) vs. 1/X plots
        i.e., Vtrans, by either interpolating the data with
        a spline function and finding X where dY/dX = 0
        that gives the most negative value of Y (opts.smooth)
        or simply the most negative value of Y (! opts.smooth)
        '''
        neg_min_x, pos_min_x = [],[]
        i = -1
        tossed = 0
        if self.opts.skipohmic:
            # Vtrans has no physical meaning for curves with negative derivatives
            self.logger.info("Skipping %s (out of %s) non-tunneling traces for Vtrans calculation." % 
                    ( len(self.ohmic), (len(self.traces))))
        
        for trace in self.avg.index.levels[0]:
            if self.opts.skipohmic and trace in self.ohmic:
                    continue
            if self.opts.nomin:
                try:
                    rootpos = self.getminroot(scipy.interpolate.UnivariateSpline(self.avg.loc[trace]['V'][self.avg.loc[trace]['V'] > 0],
                            self.avg.loc[trace]['FN'][self.avg.loc[trace]['V'] > 0], k=4, s=None))
                    rootneg = self.getminroot(scipy.interpolate.UnivariateSpline(self.avg.loc[trace]['V'][self.avg.loc[trace]['V'] < 0],
                            self.avg.loc[trace]['FN'][self.avg.loc[trace]['V'] < 0 ], k=4, s=None))
                except Exception as msg:
                    self.logger.warn("Skipped FN calculation: %s" % str(msg))
                    continue
                if np.isfinite(rootneg):
                    neg_min_x.append(rootneg)
                if np.isfinite(rootpos):
                    pos_min_x.append(rootpos)
                if np.NAN in (rootneg,rootpos):
                    self.logger.warn("No minimum found in FN derivative (-):%s, (+):%s" % (rootneg, rootpos) )
            else:
                neg_min_x.append(self.avg.loc[trace]['FN'][ self.avg[trace]['FN'][self.avg.loc[trace]['V'] < 0 ].idxmin() ])
                pos_min_x.append(self.avg.loc[trace]['FN'][ self.avg[trace]['FN'][self.avg.loc[trace]['V'] > 0 ].idxmin() ])

        if tossed: self.logger.warn("Tossed %d compliance traces during FN calculation.", tossed)
        neg_min_x = np.array(neg_min_x)
        pos_min_x = np.array(pos_min_x)
        self.FN["neg"], self.FN["pos"] = self.dohistogram(neg_min_x, "Vtrans(-)"), self.dohistogram(pos_min_x, "Vtrans(+)")

    def signedgmean(self,Y):
        '''
        Return a geometric average with the
        same sign as the input data assuming
        all values have the same sign
        '''
        if len(Y[Y<0]): Ym = -1*gmean(abs(Y))
        else: Ym = gmean(abs(Y))
        return Ym

    def gauss(self, x, *p):
        '''
        Return a gaussian function
        '''
        A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    
    def lorenz(self, x, *p):
        '''
        Return a lorenzian function
        '''
        A, mu, B = p
        return A/( (x-mu)**2 + B**2)

    def dohistogram(self, Y, label=""):
        '''
        Return a histogram of Y-values and a gaussian
        fit of the histogram, excluding values that
        exceed either the compliance limit (for current
        or current-density) or the ceiling for R. We 
        would like to include all data in the histogram,
        but outliers sometimes confuse the fitting
        routine, which defeats the purpose of machine-fitting
        '''
        
        yrange = (Y.min(),Y.max())
        if label == "J": 
            Y = Y[Y <= self.opts.compliance]
            yrange = (Y.min()-1, Y.max()+1)
        if label == "R": Y = Y[Y <= self.opts.maxr]
        if label=='DJDV': nbins = self.opts.heatmapbins
        else: nbins = self.opts.bins
        if len(Y) < 10:
            self.logger.warn("Histogram with only %d points.", len(Y))
        try:
            #TODO Why not offer density plots as an option?
            freq, bins = np.histogram(Y, range=yrange, bins=nbins, density=False)
        except ValueError as msg:
            #TODO we can now split out the file name with the bad data in it!
            self.logger.warning("Encountered this error while constructing histogram: %s", str(msg), exc_info=False)
            bins=np.array([0.,0.,0.,0.])
            freq=np.array([0.,0.,0.,0.])
        
        if len(Y):  
            Ym = self.signedgmean(Y)
            Ys = abs(Y.std())
        else:
            Ym,Ys = 0.0,0.0

        p0 = [1., Ym, Ys]
        bin_centers = (bins[:-1] + bins[1:])/2
        coeff = p0
        hist_fit = np.array([x*0 for x in range(0, len(bin_centers))])
        try:
            with self.lock:
                if self.opts.lorenzian:
                    coeff, covar = curve_fit(self.lorenz, bin_centers, freq, p0=p0, maxfev=self.opts.maxfev)
                    hist_fit = self.lorenz(bin_centers, *coeff)
                else:
                    coeff, covar = curve_fit(self.gauss, bin_centers, freq, p0=p0, maxfev=self.opts.maxfev)
                    hist_fit = self.gauss(bin_centers, *coeff)
        except RuntimeError as msg:
            if self.opts.maxfev > 100:
                self.logger.warning("|%s| Fit did not converge", label, exc_info=False)
        except ValueError as msg:
            self.logger.warning("|%s| Skipping data with ridiculous numbers in it (%s)", label, str(msg), exc_info=False )
            #coeff=p0
            #hist_fit = np.array([x*0 for x in range(0, len(bin_centers))])

        return {"bin":bin_centers, "freq":freq, "mean":coeff[1], "std":coeff[2], \
                "var":coeff[2], "bins":bins, "fit":hist_fit, "Gmean":Ym, "Gstd":Ys}

    def PrintFN(self):
        '''
        Print Vtrans values to the command line for convinience
        '''
        if self.error:
            return
        for key in ('pos', 'neg'):
            self.logger.info("|Vtrans %s| Gauss-mean: %0.4f Standard Deviation: %f" % (key, self.FN[key]['mean'], self.FN[key]['std']) )
            self.logger.info("|Vtrans %s| Geometric-mean: %0.4f Standard Deviation: %f" % (key, self.FN[key]['Gmean'], self.FN[key]['Gstd']) )