#!/usr/bin/env python3
'''
Copyright (C) 2018 Ryan Chiechi <r.c.chiechi@rug.nl>
Description:

    This is the main parsing logic for GaussFit. It was built up over
    years of trial and error with total disregard for best practices
    and is therefore a mess. However, if the warning messages are heeded,
    this software can be trusted to parse J/V data from EGaIn or CP-AFM
    reliably, as it has been tested on thousands of junctions, the results
    of have been published and they are in agreement with the litaure and
    have, to some extent, been reproduced independently.

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
# pylint: disable=line-too-long

import sys
import os
import logging
import warnings
import threading
import time
from collections import OrderedDict
from .util import signedgmean, gauss, getdistances, printFN, lorenz
from ..colors import RED, WHITE, GREEN, TEAL, YELLOW, RS
from ..logger import DelayedHandler
# import concurrent.futures

try:
    import pandas as pd
    from scipy.optimize import curve_fit, OptimizeWarning
    import scipy.interpolate
    from scipy.stats import linregress, skew,skewtest,kurtosis,kurtosistest
    import scipy.misc
    import numpy as np
    # SciPy throws a useless warning for noisy J/V traces
    warnings.filterwarnings('ignore', '.*Covariance of the parameters.*', OptimizeWarning)

except ImportError as msg:
    print("\n\t\t%s> > > Error importing numpy/pandas/scipy! %s%s%s < < <%s"
        % (RED,RS,str(msg),RED,RS))
    print("Try pip3 install <module>")
    sys.exit()

warnings.filterwarnings('ignore', '.*divide by zero.*',RuntimeWarning)
warnings.filterwarnings('ignore', '.*',UserWarning)
#warnings.filterwarnings('ignore','.*invalid value encountered in log.*',RuntimeWarning)
#warnings.filterwarnings('ignore','.*invalid value encountered in true_divide.*',RuntimeWarning)
warnings.filterwarnings('ignore', '.*invalid value encountered.*',RuntimeWarning)
warnings.filterwarnings('ignore', '.*Mean of empty slice.*', RuntimeWarning)
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

    # Class variabls
    error = False
    parsed = False
    df = pd.DataFrame()
    avg = pd.DataFrame()
    XY = OrderedDict()
    X = np.array([])
    FN = {}
    compliance_traces = []
    ohmic = []
    DJDV = {}
    GHists = {}
    NDC = []
    NDCHists = OrderedDict()
    filtered = []
    R = {}
    # traces = {}
    segments = {}
    # nofirsttrace = {}

    def __init__(self,opts,handler=None,lock=None):
        self.opts = opts
                # Pass a lock when calling inside a thread
        if lock:
            self.lock = lock
        else:
            self.lock = threading.Lock()
        # Pass your own log hanlder, e.g., when calling from a GUI
        # But make sure it supports a flush() method!
        if not handler:
            self.loghandler = DelayedHandler()
            self.loghandler.setFormatter(logging.Formatter(\
                fmt=GREEN+os.path.basename('%(name)s'+TEAL)\
                    +' %(levelname)s '+YELLOW+'%(message)s'+WHITE))
        else:
            self.loghandler = handler
        self.logger = logging.getLogger('parser')
        self.logger.addHandler(self.loghandler)
        self.logger.setLevel(getattr(logging,self.opts.loglevel.upper()))

        if not 0 < self.opts.alpha < 1:
            self.logger.error("Alpha must be between 0 and 1")
            sys.exit()

    @classmethod
    def doOutput(cls,writer):
        ''' Run through all the methods for writing output files.'''
        writer.logger.info("Writing files...")
        writer.WriteParseInfo()
        writer.WriteSummary()
        #TODO Vtrans looses its mind if you only parse positive or negative biases
        try:
            writer.WriteVtrans()
        except ValueError:
            print("Vtrans data are too broken to output ¯\_(ツ)_/¯") #pylint: disable=W1401
        writer.WriteGNUplot('Vtransplot')
        writer.WriteFN()
        writer.WriteVT()
        writer.WriteGNUplot('VTplot')
        writer.WriteGauss()
        writer.WriteSegmentedGauss()
        writer.WriteFilteredGauss()
        writer.WriteGNUplot('JVplot')
        writer.WriteGNUplot('NDCplot')
        writer.WriteData()
        writer.WriteDJDV()
        writer.WriteNDC()
        writer.WriteFiltered()
        writer.WriteData(True)
        writer.WriteLag()
        writer.WriteRData()
        writer.WriteGNUplot('Rplot')
        try:
            writer.WriteHistograms()
            writer.WriteFilteredHistograms()
            writer.WriteGNUplot('JVhistplot')
        except IndexError as msg:
            print("Error outputting histrograms %s" % str(msg))
        try:
            writer.WriteGHistogram()
        except IndexError:
            print("Error outputting Ghistrograms")
        try:
            writer.WriteGMatrix('GMatrix')
            writer.WriteGNUplot('GMatrix', ['parula.pal'])
            writer.WriteGMatrix('NDCMatrix')
            writer.WriteGNUplot('NDCMatrix')
        except IndexError:
            print("Error outputting GMatrix")
        writer.logger.info("Done!")

    def readfiles(self, fns, parse=True):
        '''Walk through input files and parse
        them into attributes '''
        frames = {}
        fns.sort()
        if isinstance(fns, str):
            fns = [fns]
        self.logger.debug('Parsing %s', ', '.join(fns))
        if self.opts.ycol > -1:
            self.logger.info("Parsing two columns of data (X=%s, Y=%s).", self.opts.xcol+1, self.opts.ycol+1)
            for f in fns:
                with open(f, 'rt') as fh:
                    _headers = fh.readline().split(self.opts.delim)
                try:
                    if _headers:
                        _x, _y = _headers[self.opts.xcol].strip(), _headers[self.opts.ycol].strip()
                        frames[f] = pd.read_csv(f,sep=self.opts.delim,
                                usecols=(_x,_y))[[_x,_y]]
                        frames[f].rename(columns = {_x:'V', _y:'J'}, inplace = True)
                        # self.logger.debug("Renaming headers %s -> V, %s -> J" % (_x, _y))
                    elif self.opts.X > self.opts.Y:
                        raise pd.errors.ParserError("xcol cannot be greater than ycol without column headers.")
                    else:
                        # self.logger.debug("No headers, manually setting V/J")
                        frames[f] = pd.read_csv(f,sep=self.opts.delim,
                                                usecols=(self.opts.xcol,
                                                        self.opts.ycol),
                                                names=('V','J'),header=0)
                except OSError as msg:
                    self.logger.warning("Skipping %s because %s", f,str(msg) )
                except pd.errors.ParserError as msg:
                    self.logger.warning("Skipping malformatted %s because %s", f,str(msg))

        else:
            self.logger.info("Parsing all columns of data.")
            for f in fns:
                try:
                    _df = pd.read_csv(f,sep=self.opts.delim,
                        index_col=self.opts.xcol,
                        header=0,
                        error_bad_lines=False,
                        warn_bad_lines=False)
                    i = 0
                    for col in _df:
                        frames['%s_%.2d' % (f,i)] = pd.DataFrame({'V':_df.index,'J':_df[col]})
                        # self.logger.debug("Adding frame %s_%.2d" % (f,i) )
                        i += 1
                except OSError as msg:
                    self.logger.warning("Skipping %s because %s" , f,str(msg))

            #self.df = pd.concat((pd.read_csv(f,sep=self.opts.delim,
            #    header=None,skiprows=1) for f in fns),ignore_index=True)
            #X,Y = [],[]
            #for row in self.df.iterrows():
            #    for y in row[1][1:]:
            #        X.append(row[1][0])
            #        Y.append(y)
            #self.df = pd.DataFrame({'V':X,'J':Y})

        if not frames:
            self.logger.error("No files to parse!")
            sys.exit()
        # Create main dataframe and parse it
        self.df = pd.concat(frames)
        # print(self.df)
        self.__parse(parse)

    def readpandas(self,df,parse):
        '''Take a pandas.DataFrame as input instead of files.'''
        self.logger.debug("Using Pandas as input")
        self.df = df
        self.__parse(parse)

    def __parse(self,parse):
        '''Read a pandas.DataFrame and compute Fowler-Nordheim
        values, log10 the J or I values and create a dictionary
        indexed by unique voltages.'''
        if (self.df.V.dtype,self.df.J.dtype) != ('float64','float64'):
            self.logger.error("Parsed data does not appear to contain numerical data!")
            self.error = True
            return
        if self.df.J.first_valid_index() is None:
            self.logger.error("Column %s is empty!", str(self.opts.ycol+1))
            self.error =True
            return
        if self.df.J.hasnans:
            self.logger.warning("Input contains non-numerical data!")
        try:
            self.df['FN'] = np.log(abs(self.df.J)/self.df.V**2)
        except ZeroDivisionError:
            self.logger.warning("Error computing FN (check your input data).")
            self.df['FN'] = np.array([x*0 for x in range(0, len(self.df['V']))])
        self.df.J.replace(0.0,value=1e-16,inplace=True)
        self.df['logJ'] = np.log10(abs(self.df.J)) # Cannot log10 zero
        self.logger.info('%s values of log|J| above compliance (%s)',
                len(self.df['logJ'][self.df['logJ']>self.opts.compliance]),self.opts.compliance)

        #The default log handler only emits when you call flush() after setDelay() called
        self.loghandler.setDelay()

        # In the event that we want to call parsing method by hand
        # we stop here when just self.df is complete
        if not parse:
            return
        self.logger.info("* * * * * * Finding segments   * * * * * * * *")
        self.loghandler.flush()
        nofirsttraces = self.__findsegments()
        self.logger.info("* * * * * * Finding traces   * * * * * * * *")
        self.loghandler.flush()
        self.__findtraces()

        if self.error:
            self.logger.error('Cannot compute statistics from these traces.')
            self.loghandler.flush()
            return # Bail if we can't parse traces

        xy = []
        for x, group in self.df.groupby('V'):
            xy.append( (x,group) )

        self.logger.info("* * * * * * Computing Lag  * * * * * * * * *")
        self.loghandler.flush()
        lag = self.dolag(xy)
        self.logger.info("* * * * * * Computing dY/dX  * * * * * * * *")
        self.loghandler.flush()
        self.__dodjdv()
        self.logger.info("* * * * * * Computing Vtrans * * * * * * * *")
        self.loghandler.flush()
        self.__findmin()
        self.logger.info("* * * * * * Computing |R|  * * * * * * * * *")
        self.loghandler.flush()
        R = self.dorect(xy)
        self.logger.info("* * * * * * Computing Gaussian  * * * * * * * * *")
        self.loghandler.flush()
        for x, group in xy:
            self.XY[x] = { 
                   "Y":group['J'],
                   "LogY":group['logJ'],
                   "hist":self.__dohistogram(group['logJ'],"J"),
                   "Y_nofirst": [0],
                   "LogY_nofirst": [0],
                   "hist_nofirst": self.__dohistogram(np.array([0]),"J"),
                   "filtered_hist":self.__dohistogram(lag[x]['filtered'],"lag"),
                   "lag":lag[x]['lagplot'],
                   "FN": group['FN'],
                   "R": R[x]
                    }
        if nofirsttraces:
            for x, group in xy:
                self.XY[x]["Y_nofirst"] = nofirsttraces[x]
                self.XY[x]["LogY_nofirst"] = [np.log10(abs(_x)) for _x in nofirsttraces[x]]
                self.XY[x]["hist_nofirst"] = self.__dohistogram(nofirsttraces[x],"J")
                    
        self.logger.info("* * * * * * Computing |V^2/J|  * * * * * * * * *")
        self.loghandler.flush()
        for x in self.XY:
            self.XY[x]['VT'] = abs(x**2/ 10**self.XY[x]['hist']['mean'] )
        self.logger.info("* * * * * * * * * * * * * * * * * * * * * * ")
        self.loghandler.flush()
        if not self.error:
            printFN(self.logger, self.FN)
        self.loghandler.unsetDelay()
        self.parsed = True

    def dorect(self,xy):
        '''
        Divide each value of Y at +V by Y at -V
        and build a histogram of rectification, R
        also construct the unique Voltage list
        '''
        r = OrderedDict()
        R = OrderedDict()
        for x, _ in xy:
            r[x] = []
        clipped = 0
        for trace in self.avg.index.levels[0]:
            for x in self.avg.loc[trace].index[self.avg.loc[trace].index >= 0]:
                if -1*x not in self.avg.loc[trace]['J']:
                    self.logger.warning("Rectification data missing voltages.")
                    if self.opts.logr:
                        r[x].append(0.)
                    else:
                        r[x].append(1.)
                    continue
                elif x == 0.0:
                    if self.opts.logr:
                        r[x].append(0.)
                    else:
                        r[x].append(1.)
                    continue
                if self.opts.logr:
                    r[x].append(np.log10(abs(self.avg.loc[trace]['J'][x]/self.avg.loc[trace]['J'][-1*x])))
                else:
                    r[x].append(abs(self.avg.loc[trace]['J'][x]/self.avg.loc[trace]['J'][-1*x]))
                if r[x][-1] > self.opts.maxr:
                    clipped += 1
        for x in reversed(list(r)):
            if x >= 0:
                R[x] = {'r':np.array(r[x]),'hist':self.__dohistogram(np.array(r[x]),"R")}
                R[-1*x] = R[x]
            if x not in R:
                self.logger.warning("Unequal +/- voltages in R-plot will be filled with R=1.")
                if self.opts.logr:
                    y = np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.])
                else:
                    y = np.array([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.])
                R[x] = {'r':y,'hist':self.__dohistogram(y,"R")}
        if clipped:
            if self.opts.logr:
                rstr = 'log|R|'
            else: rstr = '|R|'
            self.logger.info("%s values of %s exceed maxR (%s)", clipped, rstr, self.opts.maxr)
        self.logger.info("R complete.")
        return R

    def dolag(self,xy):
        '''
        Make a lag plot of Y
        '''
        lag = {}
        for x, group in xy:
            lag[x]={'lagplot':np.array([[],[]]),'filtered':np.array([])}
            Y = group['logJ']
            _lag = [[],[]]
            _filtered = []
            for i in range(0, (len(Y)-len(Y)%2), 2):
                _lag[0].append(Y[i])
                _lag[1].append(Y[i+1])
            m, b, r, _, _ = linregress(_lag[0],_lag[1])
            #self.logger.debug("R-squared: %s" % (r**2))
            distances = getdistances( (m,b), _lag[0], _lag[1] )
            min_distance = min(distances)
            self.logger.debug("Distance from lag: %s", min_distance)
            if min_distance > self.opts.lagcutoff:
                self.logger.warning("Found a high degree of scatter in lag plot (%0.4f)", min_distance)
            if r**2 < 0.9:
                self.logger.warning("Poor line-fit to lag plot (R-squared: %s)", (r**2))
            tossed = 0
            for _t in enumerate(distances):
                if _t[1] < self.opts.lagcutoff:
                    _filtered.append(_lag[0][_t[0]])
                    _filtered.append(_lag[1][_t[0]])
                else:
                    tossed += 1
            if not _filtered:
                self.logger.warning("Lag filter excluded all data at %s V", x)
            if tossed > 0:
                self.logger.info("Lag filtered excluded %s data points at %s V", tossed, x )
            lag[x]['lagplot'] = np.array(_lag)
            lag[x]['filtered'] = np.array(_filtered)
        return lag

    def __findtraces(self):
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
                    self.logger.warning(
                        'file "%s" starts and ends with weird voltages (%s -> %s)', r,s,e)
                st.append(s)
                ed.append(e)
            if self.opts.tracebyfile:
                return True
            self.logger.debug("Checking V starting from slice %s:%s",
                                traces[0][0],traces[0][1] )
            lt = len(self.df.V[ traces[0][0]:traces[0][1] ])
            for trace in traces:
                if lt != len(self.df.V[trace[0]:trace[1]]):
                    #TODO Figure out how to mask/delete these files from parsing
                    if trace[0][0] != trace[1][0]:
                        self.logger.warning(
                            'Unequal voltage steps somewhere bewteen "%s" (and) "%s"',
                                trace[0][0], trace[1][0] )
                    else:
                        self.logger.warning(
                            'Unequal voltage steps somewhere in "%s"', trace[0][0] )
                    self.loghandler.flush()
                    return False
                self.logger.debug("Trace: %s -> %s",
                        self.df.V[trace[0]],self.df.V[trace[-1]] )
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
                    # Make an iterator that aggregates elements from each of the iterables.
                    # Zip is an iterator of tuples, where the i-th tuple contains the i-th element
                    # from each of the argument sequences or iterables
                    traces.append( (t[0], t[2]) )

            except ValueError as msg:
                self.logger.warning("Did not find three-zero (EGaIn) traces!")
                self.logger.debug(str(msg))
        if not ntraces:
            self.logger.warning("This does not look like an EGaIn dataset.")
            try:
                ntraces = int(self.df.V.value_counts()[self.df.V[0]]/2) # Two end-pointss in every trace!
                for t in zip(*(iter(self.df[self.df.V == self.df.V[0]].V.index),) * 2):
                    traces.append( (t[0],t[1]) )
            except ValueError:
                self.logger.warning("Did not find three-zero (EGaIn) traces!")

        if not ntraces or not __checktraces(traces):
            # traces contains indices where traces (start, stop) in self.df.
            self.logger.warning("Recomputing traces based on repeat values")
            traces = []
            Vinit = self.df.V[0]
            trace = []
            for row in self.df[0:].iterrows():
                if row[1].V == Vinit:
                    if len(trace) <= 1:
                        trace.append(row[0])
                    elif len(trace) == 2:
                        traces.append(trace)
                        trace = [row[0]]
            ntraces = len(traces)
        if not __checktraces(traces):
            self.logger.error("Problem with traces: FN and derivative probably will not work correctly!")
            self.loghandler.flush()
        self.logger.info("Found %s traces (%s).", ntraces,len(traces) )
        idx = []
        frames = {}
        self.logger.info("Compressing forward/reverse sweeps to single traces.")
        for col,_t in enumerate(traces):
            fbtrace = self.df[traces[col][0]:traces[col][1]].sort_values('V')
            avg = OrderedDict({'J':[],'FN':[]})
            idx = []
            for x,group in fbtrace.groupby('V'):
                idx.append(x)
                avg['J'].append(signedgmean(group['J']))
                fn = np.mean(group['FN'])
                #fn = self.signedgmean(group['FN'])
                avg['FN'].append(fn)
            frames[col] = pd.DataFrame(avg,index=idx)
        try:
            self.avg = pd.concat(frames)
        except ValueError:
            self.error = True
            self.avg = pd.DataFrame()
            self.logger.error('Unable to parse traces.')
        if ntraces == 1:
            self.error = True
            self.logger.warning('Only parsed one trace!')


    def __findsegments(self):
        '''
        Break out each trace into four segments of
        0V -> Vmax, Vmax -> 0, 0V -> Vmin, Vmin -> 0V.
        '''
        #TODO set num_segments in opts
        # NOTE this is a crude hack because I forgot how Pandas works
        if self.opts.tracebyfile:
            self.logger.error("Cannot generate segments from non-EGaIn dataset.")
            return {}
        if self.opts.ycol < 0:
            self.logger.warning("Cannot find segments when all columns are parsed.")
            return {}
        try:
            if self.df.V.value_counts()[0] % 3 != 0:
                self.logger.warning("Dataset does not seem to have four segments.")
            segments = {}
            # bytrace = {}
            nofirsttrace = {}
            self.logger.info("Breaking out traces by segments of 0V -> Vmin/max.")
        except KeyError:
            self.logger.error("Could not segment data by 0's.")
            return

        for _fn in self.df.index.levels[0]:
            #TODO when parsing all columns of data, this throws:
            #  'ValueError: The truth value of a Series is ambiguous. Use a.empty, a.bool(), a.item(), a.any() or a.all().'
            # For now segmenter will just refuse to run when all columns are parsed.
            if self.df.loc[_fn]['V'][0] != 0.0:
                self.logger.warning("J/V didn't start at 0V for %s", _fn)

            _seg = 0
            _trace = 0
            _last_V = 0
            _n_traces = int(self.df.V.value_counts()[0] / 3)

            for _i in self.df.loc[_fn].index:
                J = self.df.loc[_fn]['J'][_i]
                V = self.df.loc[_fn]['V'][_i]

                if _trace > _n_traces:
                    self.logger.warning("Parsing trace %s, when there should only be %s",
                        _trace, _n_traces)

                if 0 < V > _last_V:
                    if _seg == 3:
                        _trace += 1
                    _seg = 0
                elif 0 < V <= _last_V:
                    _seg = 1
                elif 0 > V < _last_V:
                    _seg = 2
                elif 0 > V >= _last_V:
                    _seg = 3

                if _seg not in segments:
                    segments[_seg] = {'combined':{}}

                if _trace not in segments[_seg]:
                    segments[_seg][_trace] = {}

                # if _trace not in bytrace:
                #     bytrace[_trace] = {}

                if V not in segments[_seg][_trace]:
                    segments[_seg][_trace][V] = []
                
                if V not in segments[_seg]['combined']:
                    segments[_seg]['combined'][V] = []

                # if V not in bytrace[_trace]:
                #     bytrace[_trace][V] = []
                _last_V = V
                segments[_seg][_trace][V].append(J)
                segments[_seg]['combined'][V].append(J)
                # bytrace[_trace][V].append(J)

                if V not in nofirsttrace:
                    nofirsttrace[V] =[]
                if _trace > 0:
                    nofirsttrace[V].append(J)

        if len(segments.keys()) != 4:
            self.error = True
            self.logger.warning('Parsed %i segments from an EGaIn dataset!', len(segments.keys()))
        else:
            self.logger.info('Found %s segments.', len(segments.keys()))

        segmenthists = {}
        for _seg in segments:
            if _seg not in segmenthists:
                segmenthists[_seg] = {'combined':{}}
            for _trace in segments[_seg]:
                self.logger.debug('Segment: %s, Trace: %s' , _seg,_trace)
                if _trace not in segmenthists[_seg]:
                    segmenthists[_seg][_trace] = {}
                for _V in segments[_seg][_trace]:
                    if _V not in segmenthists[_seg][_trace]:
                        segmenthists[_seg][_trace][_V] = {}
                    segmenthists[_seg][_trace][_V] = self.__dohistogram(
                        np.array([np.log10(abs(_j)) for _j in segments[_seg][_trace][_V]]), label='Segmented')
            for _V in segments[_seg]['combined']:
                if _V not in segmenthists[_seg]['combined']:
                    segmenthists[_seg]['combined'][_V] = {}
                segmenthists[_seg]['combined'][_V] = self.__dohistogram(
                   np.array([np.log10(abs(_j)) for _j in segments[_seg]['combined'][_V]]), label='Segmented')

        # nofirsttracehists = {}
        for _V in nofirsttrace:
            nofirsttrace[_V] = np.array(nofirsttrace[_V])

        # tracehists = {}
        # for _trace in bytrace:
        #     if _trace not in tracehists:
        #         tracehists[_trace] = {}
        #     for _V in bytrace[_trace]:
        #         tracehists[_trace][_V] = self.__dohistogram(
        #                 np.array(bytrace[_trace][_V]), label='ByTrace')

        self.segments = segmenthists
        return nofirsttrace
        # self.bytrace = tracehists

    def __dodjdv(self):
        '''
        Fit a spline function to X/Y data,
        compute dY/dX and normalize.
        '''
        linx = np.linspace(self.df.V.min(), self.df.V.max(), 200)
        if self.opts.vcutoff > 0:
            self.logger.debug('Using %s cutoff for dj/dv', self.opts.vcutoff)
            vfilterneg,vfilterpos = np.linspace(-1*self.opts.vcutoff,0,200), np.linspace(0,self.opts.vcutoff.max(),200)
        else:
            vfilterneg,vfilterpos = np.linspace(self.df.V.min(),0,200), np.linspace(0,self.df.V.max(),200)
        if self.opts.vcutoff > 0:
            vfilterneg,vfilterpos = linx[-1*self.opts.vcutoff < linx < 0],linx[0 < linx < self.opts.vcutoff]
        else:
            vfilterneg,vfilterpos = linx[linx < 0], linx[linx > 0]

        spls = OrderedDict()
        spls_norm = OrderedDict()
        splhists = OrderedDict()
        spl_normhists = OrderedDict()
        ndc_cut, ndc_tot = 0,0
        filtered = [('Potential', 'Fit', 'Y')]
        for x in linx:
            spls[x] = []
            splhists[x] = {'spl':[],'hist':{}}
            spls_norm[x] = []
            spl_normhists[x] = {'spl':[],'hist':{}}
        for trace in self.avg.index.levels[0]:
            try:
                spl = scipy.interpolate.UnivariateSpline(
                    self.avg.loc[trace].index,self.avg.loc[trace]['J'], k=5, s=self.opts.smooth )
                dd =  scipy.interpolate.UnivariateSpline(
                    self.avg.loc[trace].index,self.avg.loc[trace]['J'], k=5, s=None).derivative(2)
            except Exception as msg: #pylint: disable=broad-except
                #TODO: Figure out how to catch all the various scipy errors here
                self.logger.error('Error in derivative calulation: %s', str(msg))
                continue
            try:
                spldd = dd(vfilterpos) #Compute d2J/dV2
                spldd += -1*dd(vfilterneg) #Compute d2J/dV2
            except ValueError as msg:
                self.logger.error('Error computing second derivative: %s', str(msg))
                continue
            if len(spldd[spldd<0]):
                    # record in the index where dY/dX is < 0 within vcutoff range
                self.ohmic.append(trace)
                if self.opts.skipohmic:
                    continue
            else:
                for row in self.avg.loc[trace].iterrows():
                    # filtered is a list containing only "clean" traces
                    filtered.append( (row[0], spl(row[0]), row[1].J) )
            err = None
            for x in spls:
                try:
                    d = spl.derivatives(x)
                except ValueError as msg:
                    err = str(msg)
                    self.logger.warning('Error computing derivative: %s', str(msg))
                    continue
                if np.isnan(d[self.opts.heatmapd]):
                    self.logger.warning("Got NaN computing dJ/dV")
                    continue
                spls[x].append(d[self.opts.heatmapd])
                splhists[x]['spl'].append(np.log10(abs(d[self.opts.heatmapd])))
                ndc = d[1] * (x/spl(x))
                spls_norm[x].append( ndc )
                ndc_tot += 1
                if 0.0 < ndc < 10.0:
                    # ndc values beyond this range can safely be considered artifacts of the numerical derivative
                    spl_normhists[x]['spl'].append( ndc )
                else:
                    # but we keep track of how many data points we toss as a sanity check
                    ndc_cut += 1
            if err:
                self.logger.error("Error while computing derivative: %s", str(err))

        self.logger.info("Non-tunneling traces: %s (out of %0d)",
                    len(self.ohmic), len(self.avg.index.levels[0]) )
        if ndc_tot:
            self.logger.info("NDC values not between 0 and 10: %s (%0.2f%%)", ndc_cut, (ndc_cut/ndc_tot)*100 )
        self.loghandler.flush()
        for x in splhists:
            splhists[x]['hist'] = self.__dohistogram(np.array(splhists[x]['spl']), label='DJDV')
            spl_normhists[x]['hist'] = self.__dohistogram(np.array(spl_normhists[x]['spl']), label='NDC')
        self.logger.info("dJdV complete.")
        self.loghandler.flush()
        self.DJDV, self.GHists, self.NDC, self.NDCHists, self.filtered = spls, splhists, spls_norm, spl_normhists, filtered

    def __findmin(self):
        '''
        Find the troughs of ln(Y/X^2) vs. 1/X plots
        i.e., Vtrans, by either interpolating the data with
        a spline function and finding X where dY/dX = 0
        that gives the most negative value of Y (opts.smooth)
        or simply the most negative value of Y (! opts.smooth)
        '''
        neg_min_x, pos_min_x = [],[]
        vmin,vmax = self.df.V.min(), self.df.V.max()
        xneg,xpos = np.linspace(vmin,0,250),np.linspace(vmax,0,250)
        tossed = 0
        if self.opts.skipohmic:
            # Vtrans has no physical meaning for curves with negative derivatives
            self.logger.info("Skipping %s (out of %s) non-tunneling traces for Vtrans calculation.",
                     len(self.ohmic), (len(self.avg.index.levels[0])) )

        for trace in self.avg.index.levels[0]:
            if self.opts.skipohmic and trace in self.ohmic:
                tossed += 1
                continue
            try:
                if not self.opts.interpolateminfn:

                    self.logger.debug('Finding FN min of plot.')
                    neg_min_x.append(self.avg.loc[trace]['FN'][ self.avg.loc[trace].index < 0 ].idxmin() )
                    pos_min_x.append(self.avg.loc[trace]['FN'][ self.avg.loc[trace].index > 0 ].idxmin() )

                else:
                    err = [False,False]
                    self.logger.debug('Finding minimum FN plot from derivative.')
                    splpos = scipy.interpolate.UnivariateSpline(np.array(self.avg.loc[trace].index[self.avg.loc[trace].index > 0]),
                                self.avg.loc[trace]['FN'][self.avg.loc[trace].index > 0].values, k=4)
                    splneg = scipy.interpolate.UnivariateSpline(np.array(self.avg.loc[trace].index[self.avg.loc[trace].index < 0]),
                                self.avg.loc[trace]['FN'][self.avg.loc[trace].index < 0].values, k=4)
                    try:
                        pos_min_x += list(np.array(splpos.derivative().roots()))
                    except ValueError as msg:
                        self.logger.warning('Error finding derivative of FN(+), falling back to linear interpolation. %s', str(msg))
                        err[0] = True
                    try:
                        neg_min_x += list(np.array(splneg.derivative().roots()))
                    except ValueError as msg:
                        self.logger.warning('Error finding derivative of FN(–), falling back to linear interpolation. %s', str(msg))
                        err[1] = True
                    if err == (False,False):
                        continue

                    self.logger.debug('Finding minimum of interpolated FN plot.')
                    splpos = scipy.interpolate.interp1d( np.array(self.avg.loc[trace].index[self.avg.loc[trace].index > 0]),
                                self.avg.loc[trace]['FN'][self.avg.loc[trace].index > 0].values,kind='linear',fill_value='extrapolate')
                    splneg = scipy.interpolate.interp1d( np.array(self.avg.loc[trace].index[self.avg.loc[trace].index < 0]),
                                self.avg.loc[trace]['FN'][self.avg.loc[trace].index < 0 ].values ,kind='linear',fill_value='extrapolate')
                    xy = {'X':[],'Y':[]}
                    for x in xneg:
                        if not np.isfinite(x):
                            continue
                        xy['Y'].append(float(splneg(x)))
                        xy['X'].append(float(x))
                    for x in xpos:
                        if not np.isfinite(x):
                            continue
                        xy['Y'].append(float(splpos(x)))
                        xy['X'].append(float(x))
                    fndf = pd.DataFrame(xy)
                    pidx = fndf[fndf.X > 0]['Y'].idxmin()
                    nidx = fndf[fndf.X < 0]['Y'].idxmin()
                    if err[0]:
                        try:
                            splpos = scipy.interpolate.UnivariateSpline(
                                fndf['X'][pidx-20:pidx+20].values,fndf['Y'][pidx-20:pidx+20].values, k=4 )
                            pos_min_x += list(np.array(splpos.derivative().roots()))
                        except Exception as msg: #pylint: disable=broad-except
                            #TODO: Figure out how to catch all the scipy errors
                            self.logger.warning('Error finding FN(+) minimum from interpolated derivative, falling back to minimum. %s', str(msg))
                            pos_min_x.append(np.mean(fndf['X'][pidx-20:pidx+20].values))
                    if err[1]:
                        try:
                            splneg = scipy.interpolate.UnivariateSpline(fndf['X'][nidx-20:nidx+20].values,fndf['Y'][nidx-20:nidx+20].values, k=4 )
                            neg_min_x += list(np.array(splneg.derivative().roots()))
                        except Exception as msg: #pylint: disable=broad-except
                            #TODO: Figure out how to catch all the scipy errors
                            self.logger.warning('Error finding FN(–) minimum from interpolated derivative, falling back to minimum. %s', str(msg))
                            neg_min_x.append(np.mean(fndf['X'][nidx-20:nidx+20].values))
            except ValueError as msg:
                self.logger.warning("Error finding minimum in trace: %s", str(msg))
        if tossed:
            self.logger.warning("Tossed %d compliance traces during FN calculation.", tossed)
        neg_min_x = np.array(list(filter(np.isfinite,neg_min_x)))
        pos_min_x = np.array(list(filter(np.isfinite,pos_min_x)))
        self.FN["neg"], self.FN["pos"] = self.__dohistogram(neg_min_x, "Vtrans(-)", True), self.__dohistogram(pos_min_x, "Vtrans(+)", True)


    def __dohistogram(self, Y, label="", density=False):
        '''
        Return a histogram of Y-values and a gaussian
        fit of the histogram, excluding values that
        exceed either the compliance limit (for current
        or current-density) or the ceiling for R. We
        would like to include all data in the histogram,
        but outliers sometimes confuse the fitting
        routine, which defeats the purpose of machine-fitting
        '''

        try:
            yrange = (Y.min(),Y.max())
        except ValueError as msg:
            self.logger.error("Error ranging data for histogram: %s", str(msg))
            yrange = (0,0)

        if label == "J" or label == "lag":
            Y = Y[Y <= self.opts.compliance]
            if yrange != (0,0):
                yrange = (Y.min()-1, Y.max()+1)
        if label == "R":
            Y = Y[Y <= self.opts.maxr]
        if label in ('DJDV','NDC'):
            nbins = self.opts.heatmapbins
        else: nbins = self.opts.bins
        if len(Y) < 10:
            self.logger.warning("Histogram with only %d points.", len(Y))
        try:
            #TODO Why not offer density plots as an option?
            freq, bins = np.histogram(Y, range=yrange, bins=nbins, density=density)
        except ValueError as msg:
            #TODO we can now split out the file name with the bad data in it!
            self.logger.warning("Encountered this error while constructing histogram: %s", str(msg), exc_info=False)
            bins=np.array([0.,0.,0.,0.])
            freq=np.array([0.,0.,0.,0.])

        if len(Y):
            Ym = signedgmean(Y)
            Ys = abs(Y.std())
        else:
            Ym,Ys = 0.0,0.0

        p0 = [1., Ym, Ys]
        bin_centers = (bins[:-1] + bins[1:])/2
        coeff = p0
        covar = None #pylint: disable=unused-variable
        hist_fit = np.array([x*0 for x in range(0, len(bin_centers))])
        try:
            with self.lock:
                if self.opts.lorenzian:
                    coeff, covar = curve_fit(lorenz, bin_centers, freq, p0=p0, maxfev=self.opts.maxfev) #pylint: disable=unbalanced-tuple-unpacking
                    hist_fit = lorenz(bin_centers, *coeff)
                else:
                    coeff, covar = curve_fit(gauss, bin_centers, freq, p0=p0, maxfev=self.opts.maxfev) #pylint: disable=unbalanced-tuple-unpacking
                    hist_fit = gauss(bin_centers, *coeff)
        except RuntimeError as msg:
            if self.opts.maxfev > 100:
                self.logger.warning("|%s| Fit did not converge", label, exc_info=False)
        except ValueError as msg:
            self.logger.warning("|%s| Skipping data with ridiculous numbers in it (%s)", label, str(msg), exc_info=False )
            #coeff=p0

        # if covar is not None:
        #     try:
        #         self.logger.debug('Covariance: %s ', np.sqrt(np.diag(covar)) )
        #     except FloatingPointError:
        #         self.logger.debug('Could not compute sqrt of covariance: %s', np.diag(covar))

        skewstat, skewpval = skewtest(freq)
        kurtstat, kurtpval = kurtosistest(freq)
        return {"bin":bin_centers, "freq":freq, "mean":coeff[1], "std":coeff[2], \
                "var":coeff[2], "bins":bins, "fit":hist_fit, "Gmean":Ym, "Gstd":Ys,\
                "skew":skew(freq), "kurtosis":kurtosis(freq), "skewstat":skewstat, "skewpval":skewpval,
                "kurtstat":kurtstat, "kurtpval":kurtpval}

    def wait(self):
        '''
        Wait at most 60 seconds for either an error to occur or
        for the parser to complete.
        '''
        self.logger.debug("Waiting for parser to complete.")
        t = 0
        while not self.parsed and not self.error:
            if t > 60:
                self.logger.error("Timeout waiting for parser to complete.")
                sys.exit(-1)
                break
            time.sleep(1)
            t += 1
        if self.error:
            self.logger.error("!!! Parser completing with error, check the results carefully !!!")
            self.loghandler.flush()

    def getXY(self):
        '''Public method to return the main XY dictionary'''
        self.wait()
        if self.error:
            return {}
        else:
            return self.XY

    def getFN(self):
        '''Public method to reutrn the main FN dictionary'''
        self.wait()
        if self.error:
            return {}
        else:
            return self.FN
