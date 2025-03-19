#!/usr/bin/env python3
'''
Copyright (C) 2023 Ryan Chiechi <ryan.chiechi@ncsu.edu>
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

import asyncio
import sys
import os
import logging
from logging.handlers import QueueListener, QueueHandler
from multiprocessing import Queue
import warnings
import time
import pickle
from collections import OrderedDict
from gaussfit.args import VERSION
from gaussfit.parse.libparse.util import printFN, throwimportwarning
from SLM.extract import findvtrans, findG
from SLM.util import SLM_func, Gamma_func, slm_param_func
from gaussfit.colors import WHITE, GREEN, TEAL, YELLOW
from gaussfit.logger import DelayedHandler
from gaussfit.parse.libparse.util import gettmpfilename
from gaussfit.parse.libparse.util import getfilechecksum
from gaussfit.parse.libparse.dohistogram import dohistogram
from gaussfit.parse.libparse.findsegments import findsegments
from gaussfit.parse.libparse.dolag  import dolag
import platform
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
warnings.filterwarnings('error', '.*Degrees of freedom <= 0 for slice.*', RuntimeWarning)
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
    from gaussfit.parse.libparse import dodjdv
    from gaussfit.parse.libparse import findmin, old_findmin
    from gaussfit.parse.libparse import dorect
    from gaussfit.parse.libparse import doconductance
    from gaussfit.parse.libparse import doslm

    # Class variables
    VERSION = '1.0.2a'
    error = False
    parsed = False
    file_hashes = {}
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
    G = {}  # Conductance indexed by trace
    SLM = {'G': {},
           'Vtpos': {},
           'Vtneg': {},
           'epsillon': {},
           'gamma': {},
           'big_gamma': {},
           'fit': {},
           'calc': {},
           'exp': {},
           'full': {}}  # SLM inputs and outputs by trace
    segments = {}
    segments_nofirst = {}
    logger = logging.getLogger(__package__)
    logqueue = Queue()

    def __init__(self, opts, handler=None):
        self.opts = opts
        # Pass your own log handler, e.g., when calling from a GUI
        # But make sure it supports flush(), setDelay() and unsetDelay() methods!
        if not handler:
            self.loghandler = DelayedHandler()
            self.loghandler.setFormatter(logging.Formatter(
                fmt=GREEN + os.path.basename(
                    '%(name)s' + TEAL) + ' %(levelname)s ' + YELLOW + '%(message)s' + WHITE))
        else:
            self.loghandler = handler
        self.logger.addHandler(self.loghandler)
        self.logger.setLevel(getattr(logging, self.opts.loglevel.upper()))
        self.loglistener = QueueListener(self.logqueue, self.loghandler)
        self.loglistener.start()
        self.logger.info("Gaussfit v%s", VERSION)

        if not 0 < self.opts.alpha < 1:
            self.logger.error("Alpha must be between 0 and 1")
            sys.exit()

    def _checkfordupe(self, f):
        _dupe = False
        _digest = getfilechecksum(f)
        if _digest in self.file_hashes:
            self.logger.warning(f'{self.file_hashes[_digest]} and {f} are identical!')
            _dupe = True
        self.file_hashes[_digest] = f
        return _dupe

    def _dedupefiles(self, _fns):
        fns = []
        if isinstance(_fns, str):
            _fns = [_fns]
        for _f in _fns:
            if self._checkfordupe(_f):
                self.logger.warning(f'Refusing to parse duplicate file {_f}.')
                self.logger.warning('Parsing identical files will skew the results.')
            else:
                fns.append(_f)
        fns.sort()
        return fns

    async def readfiles(self, _fns, parse=True):
        '''Walk through input files and parse
        them into attributes '''

        frames = {}
        fns = self._dedupefiles(_fns)

        self.logger.debug('Parsing %s', ', '.join(fns))
        if self.opts.ycol > -1:
            self.logger.info("Parsing two columns of data (X=%s, Y=%s).", self.opts.xcol + 1, self.opts.ycol + 1)
            for f in fns:
                with open(f, 'rb') as fh:
                    try:
                        _headers = fh.readline().split(bytes(self.opts.delim, encoding=self.opts.encoding))
                        _headers = list(map(lambda x: str(x, encoding=self.opts.encoding), _headers))
                    except UnicodeDecodeError:
                        self.logger.warning("Encountered an illegal unicode character in headers.")
                try:
                    if _headers:
                        _x, _y = _headers[self.opts.xcol].strip(), _headers[self.opts.ycol].strip()
                        frames[f] = pd.read_csv(f, sep=self.opts.delim, encoding=self.opts.encoding,
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

        if not frames:
            self.logger.error("No files to parse!")
            sys.exit()
        # Create main dataframe and parse it
        self.df = pd.concat(frames)
        if self.opts.xrange > 0:
            self.logger.info(f"Pruning x-axis to +/- {self.opts.xrange}.")
            self.df = self.df[self.df.V > -1 * abs(self.opts.xrange)]
            self.df = self.df[self.df.V < abs(self.opts.xrange)]
        # print(self.df)
        await self.__parse(parse)

    def readpandas(self, df, parse):
        '''Take a pandas.DataFrame as input instead of files.'''
        self.logger.debug("Using Pandas as input")
        self.df = df
        self.__parse(parse)

    async def __parse(self, parse):
        '''Read a pandas.DataFrame and compute Fowler-Nordheim
        values, log10 the J or I values and create a dictionary
        indexed by unique voltages.'''

        if (self.df.V.dtype, self.df.J.dtype) != ('float64', 'float64'):
            self.logger.error("Parsed data does not appear to contain numerical data!")
            self.error = True
            return
        if self.df.J.first_valid_index() is None:
            self.logger.error("Column %s is empty!", str(self.opts.ycol + 1))
            self.error = True
            return
        if self.df.J.hasnans:
            self.logger.warning("Input contains non-numerical data!")
        try:
            self.df['FN'] = np.log(abs(self.df.J) / self.df.V ** 2)
        except ZeroDivisionError:
            self.logger.warning("Error computing FN (check your input data).")
            self.df['FN'] = np.array([x * 0 for x in range(0, len(self.df['V']))])
        self.df.replace({'J':0.0}, value=1e-16, inplace=True)
        self.df['logJ'] = np.log10(abs(self.df.J))  # Cannot log10 zero
        self.logger.info('%s values of log|J| above compliance (%s)',
                         len(self.df['logJ'][self.df['logJ'] > self.opts.compliance]), self.opts.compliance)

        # The default log handler only emits when you call flush() after setDelay() called
        self.loghandler.setDelay()

        # In the event that we want to call parsing method by hand
        # we stop here when just self.df is complete
        if parse:
            await self.__parsedataset()

    async def __parsedataset(self):
        children = []
        tasks = []
        xy = []
        for x, group in self.df.groupby('V'):
            xy.append((x, group))
        self.logger.info("* * * * * * Finding segments   * * * * * * * *")
        tasks.append(asyncio.create_task(findsegments(self.df, que=self.logqueue), name="findsegments"))
        self.logger.info("* * * * * * Finding traces   * * * * * * * *")
        self.loghandler.flush()
        self.findtraces()
        self.SLM['G_avg'] = self.doconductance()
        self.logger.info("* * * * * * Computing Lag  * * * * * * * * *")
        tasks.append(asyncio.create_task(dolag(xy, que=self.logqueue), name="dolag"))
        self.dodjdv()
        if self.opts.oldfn:
            self.old_findmin()
        else:
            self.findmin()
        self.SLM['Vtposavg'] = self.FN["pos"]
        self.SLM['Vtnegavg'] = self.FN["neg"]
        R = self.dorect(xy)
        results = await asyncio.gather(*tasks, return_exceptions=True)
        # Map results to task names
        task_results = {
            task.get_name(): result 
            for task, result in zip(tasks, results)
        }
        self.error, self.segments, self.segmenthists_nofirst, nofirsttrace = task_results["findsegments"]
        lag = task_results["dolag"]
        self.logger.info("* * * * * * Computing Gaussians  * * * * * * * * *")
        self.loghandler.flush()
        for x, group in xy:
            self.XY[x] = {
                "Y": group['J'],
                "LogY": group['logJ'],
                "hist": dohistogram(group['logJ'], label="J", warnings=True, que=self.logqueue),
                "Y_nofirst": [0],
                "LogY_nofirst": [0],
                "hist_nofirst": dohistogram(np.array([0]), label="J", que=self.logqueue),
                "filtered_hist": dohistogram(lag[x]['filtered'], label="lag", que=self.logqueue),
                "lag": lag[x]['lagplot'],
                "FN": group['FN'],
                "R": R[x]}
        if self.opts.heatmapd == 0:
            self.GHists = OrderedDict()
            for x in self.XY:
                self.GHists[x] = {}
                self.GHists[x]['hist'] = self.XY[x]['hist']
        if nofirsttrace:
            for x, group in xy:
                self.XY[x]["Y_nofirst"] = nofirsttrace[x]
                self.XY[x]["LogY_nofirst"] = [np.log10(abs(_x)) for _x in nofirsttrace[x]]
                self.XY[x]["hist_nofirst"] = dohistogram(nofirsttrace[x], label="J", que=self.logqueue)

        self.logger.info("* * * * * * Computing |V^2/J|  * * * * * * * * *")
        self.loghandler.flush()
        for x in self.XY:
            self.XY[x]['VT'] = abs(x**2 / 10**self.XY[x]['hist']['mean'])
        self.logger.info("* * * * * * Computing SLM from Gaussian LogJ  * * * * * * * * *")
        self.loghandler.flush()
        _v, _j = [], []
        for x in self.XY:
            _v.append(x)
            _j.append(self.XY[x]['hist']['mean'])
        SLM_logger = logging.getLogger(__package__+".SLM")
        SLM_logger.addHandler(QueueHandler(self.logqueue))
        self.SLM['Gauss'] = {}
        self.SLM['Gauss']['FN'] = findvtrans(_v, _j, logger=SLM_logger, unlog=True)
        self.SLM['Gauss']['G'] = findG(_v, _j, logger=SLM_logger, unlog=True)['slope']
        epsillon, gamma = slm_param_func(self.SLM['Gauss']['FN']['vt_pos'], self.SLM['Gauss']['FN']['vt_neg'])
        big_gamma = Gamma_func(self.SLM['Gauss']['G'], self.opts.nmolecules,
                               self.SLM['Gauss']['FN']['vt_pos'], self.SLM['Gauss']['FN']['vt_neg'])
        self.SLM['Gauss']['epsillon'] = epsillon
        self.SLM['Gauss']['gamma'] = gamma
        self.SLM['Gauss']['big_gamma'] = big_gamma
        self.logger.info("* * * * * * Computing SLM  * * * * * * * * *")
        self.loghandler.flush()
        self.logger.info(f"Fit {self.doslm()} traces to SLM.")
        self.logger.info("* * * * * * * * * * * * * * * * * * * * * * ")
        self.loghandler.flush()
        if not self.error:
            printFN(self.logger, self.FN)
        self.logger.info("Vtrans +/- from Gaussian LogJ data:")
        self.logger.info(f"{self.SLM['Gauss']['FN']['vt_pos']:0.2f} / {self.SLM['Gauss']['FN']['vt_neg']:0.2f}")
        self.logger.info("SLM from Gaussian LogJ data:")
        self.logger.info(f"G = {self.SLM['Gauss']['G']:0.2E}, ε = {epsillon:0.2f}, γ = {gamma:0.2f}, Γ = {big_gamma:0.2E}")
        self.loghandler.unsetDelay()

        try:
            self.loglistener.stop()
        except EOFError:
            pass

        if self.error:
            self.logger.error('Cannot compute statistics from these traces. (Did you set segments correctly?)')
            if self.opts.force:
                self.logger.warn('Continuing anyway.')
                self.error = False
            self.loghandler.flush()
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
