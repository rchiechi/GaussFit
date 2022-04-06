#!/usr/bin/env python3
'''
Copyright (C) 2022 Ryan Chiechi <ryan.chiechi@ncsu.edu>
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
from logging.handlers import QueueListener, QueueHandler
import warnings
# import threading
import time
import pickle
from tempfile import NamedTemporaryFile
from collections import OrderedDict
from gaussfit.parse.libparse.util import printFN, throwimportwarning
from gaussfit.colors import WHITE, GREEN, TEAL, YELLOW
from gaussfit.logger import DelayedHandler
from gaussfit.parse.libparse.dummies import dummyListener, dummyPopen
# import concurrent.futures
from multiprocessing import Process, Queue
import platform  # avoids TypeError: cannot pickle '_thread.lock' object error
if platform.system() == "Darwin":
    from multiprocessing import set_start_method
    set_start_method("fork")
elif platform.system() == "Windows":
    from multiprocessing import set_start_method
    set_start_method("spawn")

try:
    import pandas as pd
    from scipy.optimize import OptimizeWarning
    import numpy as np
    # SciPy throws a useless warning for noisy J/V traces
    warnings.filterwarnings('ignore', '.*Covariance of the parameters.*', OptimizeWarning)

except ImportError as msg:
    throwimportwarning(msg)

warnings.filterwarnings('ignore', '.*divide by zero.*', RuntimeWarning)
warnings.filterwarnings('ignore', '.*', UserWarning)
# warnings.filterwarnings('ignore','.*invalid value encountered in log.*',RuntimeWarning)
# warnings.filterwarnings('ignore','.*invalid value encountered in true_divide.*',RuntimeWarning)
warnings.filterwarnings('ignore', '.*invalid value encountered.*', RuntimeWarning)
warnings.filterwarnings('ignore', '.*Mean of empty slice.*', RuntimeWarning)
# warnings.filterwarnings('error','.*Degrees of freedom <= 0 for slice.*', RuntimeWarning)
# warnings.filterwarnings('ignore','.*impossible result.*',UserWarning)


class Parse():
    '''
    This is the main parsing class that takes input data
    in the form of text files, parses them into dictionary
    and list attributes and provides methods for performing
    operations: Gaussian fits, F-N calculations, Rectification,
    Vtrans, and dJ/DV.
    '''

    from gaussfit.parse.libparse import findtraces
    from gaussfit.parse.libparse import findsegments
    from gaussfit.parse.libparse import dodjdv
    from gaussfit.parse.libparse import findmin
    from gaussfit.parse.libparse import dorect
    from gaussfit.parse.libparse import dohistogram
    from gaussfit.parse.libparse import dolag

    # Class variables
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
    GHists = OrderedDict()
    NDC = []
    NDCHists = OrderedDict()
    filtered = []
    R = {}
    segments = {}
    segments_nofirst = {}
    # called_from_gui = False
    loglistener = None
    logger = logging.getLogger('parser')

    def __init__(self, opts, handler=None):
        self.opts = opts
        # Pass your own log hanlder, e.g., when calling from a GUI
        # But make sure it supports flush(), setDelay() and unsetDelay() methods!
        if not handler:
            self.loghandler = DelayedHandler()
            self.loghandler.setFormatter(logging.Formatter(
                fmt=GREEN+os.path.basename(
                    '%(name)s'+TEAL)+' %(levelname)s '+YELLOW+'%(message)s'+WHITE))
            logqueue = Queue(-1)
            queuehanlder = QueueHandler(logqueue)
            self.loglistener = QueueListener(logqueue, self.loghandler)
            self.logger.addHandler(queuehanlder)
            self.loglistener.start()
        else:
            self.loghandler = handler
            self.loglistener = dummyListener()
            self.logger.addHandler(self.loghandler)
        self.logger.setLevel(getattr(logging, self.opts.loglevel.upper()))

        if not 0 < self.opts.alpha < 1:
            self.logger.error("Alpha must be between 0 and 1")
            sys.exit()

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
                        frames[f] = pd.read_csv(f, sep=self.opts.delim,
                                                usecols=(_x, _y))[[_x, _y]]
                        frames[f].rename(columns={_x: 'V', _y: 'J'}, inplace=True)
                        # self.logger.debug("Renaming headers %s -> V, %s -> J" % (_x, _y))
                    elif self.opts.X > self.opts.Y:
                        raise pd.errors.ParserError("xcol cannot be greater than ycol without column headers.")
                    else:
                        # self.logger.debug("No headers, manually setting V/J")
                        frames[f] = pd.read_csv(f, sep=self.opts.delim,
                                                usecols=(self.opts.xcol,
                                                         self.opts.ycol),
                                                names=('V', 'J'), header=0)
                except OSError as msg:
                    self.logger.warning("Skipping %s because %s", f, str(msg))
                except pd.errors.ParserError as msg:
                    self.logger.warning("Skipping malformatted %s because %s", f, str(msg))

        else:
            self.logger.info("Parsing all columns of data.")
            for f in fns:
                try:
                    _df = pd.read_csv(f, sep=self.opts.delim,
                                      index_col=self.opts.xcol,
                                      header=0,
                                      error_bad_lines=False,
                                      warn_bad_lines=False)
                    i = 0
                    for col in _df:
                        frames['%s_%.2d' % (f, i)] = pd.DataFrame({'V': _df.index, 'J': _df[col]})
                        # self.logger.debug("Adding frame %s_%.2d" % (f,i) )
                        i += 1
                except OSError as msg:
                    self.logger.warning("Skipping %s because %s", f, str(msg))

            # self.df = pd.concat((pd.read_csv(f,sep=self.opts.delim,
            #    header=None,skiprows=1) for f in fns),ignore_index=True)
            # X,Y = [],[]
            # for row in self.df.iterrows():
            #    for y in row[1][1:]:
            #        X.append(row[1][0])
            #        Y.append(y)
            # self.df = pd.DataFrame({'V':X,'J':Y})

        if not frames:
            self.logger.error("No files to parse!")
            sys.exit()
        # Create main dataframe and parse it
        self.df = pd.concat(frames)
        # print(self.df)
        self.__parse(parse)

    def readpandas(self, df, parse):
        '''Take a pandas.DataFrame as input instead of files.'''
        self.logger.debug("Using Pandas as input")
        self.df = df
        self.__parse(parse)

    def __parse(self, parse):
        '''Read a pandas.DataFrame and compute Fowler-Nordheim
        values, log10 the J or I values and create a dictionary
        indexed by unique voltages.'''

        if (self.df.V.dtype, self.df.J.dtype) != ('float64', 'float64'):
            self.logger.error("Parsed data does not appear to contain numerical data!")
            self.error = True
            return
        if self.df.J.first_valid_index() is None:
            self.logger.error("Column %s is empty!", str(self.opts.ycol+1))
            self.error = True
            return
        if self.df.J.hasnans:
            self.logger.warning("Input contains non-numerical data!")
        try:
            self.df['FN'] = np.log(abs(self.df.J)/self.df.V**2)
        except ZeroDivisionError:
            self.logger.warning("Error computing FN (check your input data).")
            self.df['FN'] = np.array([x*0 for x in range(0, len(self.df['V']))])
        self.df.J.replace(0.0, value=1e-16, inplace=True)
        self.df['logJ'] = np.log10(abs(self.df.J))  # Cannot log10 zero
        self.logger.info('%s values of log|J| above compliance (%s)',
                         len(self.df['logJ'][self.df['logJ'] > self.opts.compliance]), self.opts.compliance)

        # The default log handler only emits when you call flush() after setDelay() called
        self.loghandler.setDelay()

        # In the event that we want to call parsing method by hand
        # we stop here when just self.df is complete
        if parse:
            self.__parsedataset()

    def __parsedataset(self):
        children = []
        # if self.called_from_gui:
        parent_conn = NamedTemporaryFile()  # multiprocess Pipe is not thread safe
        child_conn = parent_conn
        # else:
        #    parent_conn, child_conn = Pipe(duplex=False)
        if platform.system() == 'Linux':
            __p = Process(target=self.findsegments, args=(child_conn,))
            __p.start()
            children.append([parent_conn, __p])
        else:
            self.findsegments(child_conn)
            children.append([parent_conn, dummyPopen()])
        self.logger.info("* * * * * * Finding traces   * * * * * * * *")
        self.loghandler.flush()
        self.findtraces()

        xy = []
        for x, group in self.df.groupby('V'):
            xy.append((x, group))

        # if self.called_from_gui:
        parent_conn = NamedTemporaryFile()
        child_conn = parent_conn
        # else:
        #    parent_conn, child_conn = Pipe(duplex=False)
        if platform.system() == 'Linux':
            __p = Process(target=self.dolag, args=(child_conn, xy,))
            __p.start()
            children.append([parent_conn, __p])
        else:
            self.dolag(child_conn, xy)
            children.append([parent_conn, dummyPopen()])
        self.dodjdv()
        self.findmin()
        R = self.dorect(xy)

        #if self.called_from_gui:
        children[1][1].join()
        lag = pickle.load(children[1][0])
        #else:
        #    lag = children[1][0].recv()
        #    children[1][1].join()
        self.logger.info("* * * * * * Computing Gaussians  * * * * * * * * *")
        self.loghandler.flush()
        for x, group in xy:
            self.XY[x] = {
                "Y": group['J'],
                "LogY": group['logJ'],
                "hist": self.dohistogram(group['logJ'], "J"),
                "Y_nofirst": [0],
                "LogY_nofirst": [0],
                "hist_nofirst": self.dohistogram(np.array([0]), "J"),
                "filtered_hist": self.dohistogram(lag[x]['filtered'], "lag"),
                "lag": lag[x]['lagplot'],
                "FN": group['FN'],
                "R": R[x]}
        if self.opts.heatmapd == 0:
            self.GHists = OrderedDict()
            for x in self.XY:
                self.GHists[x] = {}
                self.GHists[x]['hist'] = self.XY[x]['hist']
        # if self.called_from_gui:
        children[0][1].join()
        self.error, self.segments, self.segmenthists_nofirst, nofirsttrace = pickle.load(children[0][0])
        # else:
        #    self.error, self.segments, self.segmenthists_nofirst, nofirsttrace = children[0][0].recv()
        #    children[0][1].join()
        if nofirsttrace:
            for x, group in xy:
                self.XY[x]["Y_nofirst"] = nofirsttrace[x]
                self.XY[x]["LogY_nofirst"] = [np.log10(abs(_x)) for _x in nofirsttrace[x]]
                self.XY[x]["hist_nofirst"] = self.dohistogram(nofirsttrace[x], "J")

        self.logger.info("* * * * * * Computing |V^2/J|  * * * * * * * * *")
        self.loghandler.flush()
        for x in self.XY:
            self.XY[x]['VT'] = abs(x**2 / 10**self.XY[x]['hist']['mean'])
        self.logger.info("* * * * * * * * * * * * * * * * * * * * * * ")
        self.loghandler.flush()
        if not self.error:
            printFN(self.logger, self.FN)
        self.loghandler.unsetDelay()

        __listener = getattr(self, 'loglistener')
        if callable(__listener):
            self.loglistener.stop()

        if self.error:
            self.logger.error('Cannot compute statistics from these traces. (Did you set segments correctly?)')
            self.loghandler.flush()
            return  # Bail if we can't parse traces
        self.parsed = True

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
