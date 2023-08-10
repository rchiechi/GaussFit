'''
Copyright (C) 2022 Ryan Chiechi <ryan.chiechi@ncsu.edu>
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
import datetime
import csv
import warnings
from shutil import copyfile
import logging
from gaussfit.colors import GREEN, TEAL, YELLOW, WHITE
from gaussfit.args import Opts as opts

logger = logging.getLogger('output')
loghandler = logging.StreamHandler()
loghandler.setFormatter(logging.Formatter(
    fmt=GREEN + os.path.basename('%(name)s' + TEAL) + ' %(levelname)s ' + YELLOW + '%(message)s' + WHITE))
logger.addHandler(loghandler)

warnings.filterwarnings('ignore', '.*comparison.*', FutureWarning)


class Writer():
    '''The main Writer class for creating text files of parsed data.'''
    from gaussfit.output.libwriter import WriteHistograms
    from gaussfit.output.libwriter import WriteFilteredHistograms
    from gaussfit.output.libwriter import WriteSegmentedGauss
    from gaussfit.output.libwriter import WriteVtrans
    from gaussfit.output.libwriter import WriteFN
    from gaussfit.output.libwriter import WriteFilteredGauss
    from gaussfit.output.libwriter import WriteGauss
    from gaussfit.output.libwriter import WriteVT
    from gaussfit.output.libwriter import WriteSLM
    from gaussfit.output.libwriter import WriteData
    from gaussfit.output.libwriter import WriteDJDV
    from gaussfit.output.libwriter import WriteNDC
    from gaussfit.output.libwriter import WriteFiltered
    from gaussfit.output.libwriter import WriteLag
    from gaussfit.output.libwriter import WriteRData
    from gaussfit.output.libwriter import WriteGHistogram
    from gaussfit.output.libwriter import WriteGMatrix

    def __init__(self, parser):
        self.parser = parser
        self.opts = opts
        if not os.path.exists(opts.out_dir):
            logger.info("Creating %s", opts.out_dir)
            os.makedirs(opts.out_dir)
        if not os.path.exists(opts.slm_dir):
            logger.info("Creating %s", opts.slm_dir)
            os.makedirs(opts.slm_dir)

    def __getattr__(self, name):
        try:
            return getattr(self.parser, name)  # inheret the methods of self.parser
        except AttributeError as msg:
            raise AttributeError("Writer object has no attribute '%s'" % name) from msg

    def WriteParseInfo(self, extra=''):
        '''Write some summary information about the parameters
        used to parse the input data.'''
        _fn = os.path.join(self.opts.out_dir, self.opts.outfile + "_parseinfo.txt")
        with open(_fn, 'a') as _fh:
            _fh.write("Parsed: %s\n***\n" % str(datetime.datetime.today().ctime()))
            _t = str(vars(self.opts))
            _t = _t.replace(",", "\n").replace("[", "\n[")
            _fh.write(_t)
            _fh.write(extra + "\n")
            _fh.write("\n***\n")

    def WriteSummary(self):
        '''Write a summary of the traces that were parsed.'''
        try:
            _fn = os.path.join(self.opts.out_dir, self.opts.outfile + "_Summary.txt")
            self.df.to_csv(_fn, sep=self.opts.delim)
        except AttributeError:
            logger.warning("No derivative data to summarize")
        try:
            _fn = os.path.join(self.opts.out_dir, self.opts.outfile + "_Traces.txt")
            self.avg.to_csv(_fn, sep=self.opts.delim)
        except AttributeError:
            logger.warning("No averaged data to summarize")

    def WriteGNUplot(self, gpinbn, tocopy=None):
        '''Write GNUPlot sciprts for plotting the ASCII files output by the other writers.'''
        tocopy = tocopy or []
        absdir = os.path.dirname(os.path.abspath(__file__))
        tdir = os.path.join(absdir, '../../templates/')
        gpintp = os.path.join(tdir, gpinbn + '.gpin')
        try:
            nsub = str(open(gpintp, 'rt').read()).count('%s')
        except FileNotFoundError:
            logger.warning("Could not read template file %s", gpintp)
            return
        ssub = []
        for _ in range(nsub):
            ssub.append(self.opts.outfile)
        txt = open(gpintp, 'rt').read() % tuple(ssub)
        _fh = open(os.path.join(self.opts.out_dir, f'{self.opts.outfile}_{gpinbn}.gpin'), 'wt')
        _fh.write(txt)
        _fh.close()
        for _fn in tocopy:
            copyfile(os.path.join(tdir, _fn),
                     os.path.join(self.opts.out_dir, _fn))

    def WriteGeneric(self, dataset, bfn, labels=None):
        '''Write a generic set of data expecting an n-dimensional array'''
        labels = labels or []
        if labels and len(labels) != len(dataset):
            logger.error("Length of column labels does not match number of data columns for WriteGeneric!")
            return

        lencola = len(dataset[0])
        for d in dataset:
            if len(d) != lencola:
                logger.error("Length of columns differs for WriteGeneric!")
                return

        _fn = os.path.join(self.opts.out_dir, f"{self.opts.outfile}_{bfn}.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            if len(labels):
                writer.writerow(labels)

            for n in range(0, lencola):
                row = []
                for _t in enumerate(dataset):
                    row.append(_t[1][n])
                writer.writerow(row)
