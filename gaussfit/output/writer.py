#!/usr/bin/env python3
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
#pylint: disable=line-too-long

import os
import warnings
import datetime
import csv
from shutil import copyfile
import logging
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
except ImportError as msg:
    pass #We catch numpy import errors in Parser.py
warnings.filterwarnings('ignore','.*comparison.*',FutureWarning)



class Writer():
    '''The main Writer class for creating text files of parsed data.'''
    def __init__(self,parser):
        self.parser = parser
        self.opts = self.parser.opts
        if not os.path.exists(parser.opts.out_dir):
            logger.info("Creating %s" , parser.opts.out_dir)
            os.mkdir(parser.opts.out_dir)

        # self.df = len(self.opts.in_files)-1 or 1

    def __getattr__(self, name):
        try:
            return getattr(self.parser, name) # 'inheret' the methods of self.parser
        except AttributeError as msg:
            raise AttributeError("Writer object has no attribute '%s'" % name) from msg

    def WriteParseInfo(self,extra=''):
        '''Write some summary information about the parameters
        used to parse the input data.'''
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_parseinfo.txt")
        with open(_fn, 'a') as _fh:
            _fh.write("Parsed: %s\n***\n" % str(datetime.datetime.today().ctime()) )
            _t = str(vars(self.opts))
            _t = _t.replace(",","\n").replace("[","\n[")
            _fh.write(_t)
            _fh.write(extra+"\n")
            _fh.write("\n***\n")

    def WriteSummary(self):
        '''Write a summary of the traces that were parsed.'''
        try:
            _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Summary.txt")
            self.df.to_csv(_fn,sep=self.opts.delim)
        except AttributeError:
            logger.warning("No derivative data to summarize")
        try:
            _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Traces.txt")
            self.avg.to_csv(_fn,sep=self.opts.delim)
        except AttributeError:
            logger.warning("No averaged data to summarize")

    def WriteHistograms(self):
        '''Write all of the underlying histograms and associated statistics used to compute
           Gaussian mean and variance of J, R and Vtrans from the raw input data.'''
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Histograms.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = []
            #for x in self.XY: headers += ["Log |J| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x, \
            #        "Skew (%0.4f)"%x, "Kurtosis (%0.4f)"%x, "Skew test (%0.4f)"%x, "Skew pvalue (%0.4f)"%x]
            for x in self.XY:
                headers += ["Log |J| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x]
            writer.writerow(headers)
            for i in range(0, len( self.XY[list(self.XY.keys())[0]]['hist']['bin'] ) ):
                row = []
                for x in self.XY:
                    row += ["%0.4f"%self.XY[x]['hist']['bin'][i],
                        "%s"%self.XY[x]['hist']['freq'][i],
                        "%0.4f"%self.XY[x]['hist']['fit'][i]]
                writer.writerow(row)

        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Histograms_stats.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = ["Voltage", "Skew", "Kurtosis", "Skew zscore", "Skew pvalue", "Kurtosis zscore", "Kurtosis pvalue"]
            writer.writerow(headers)
            for x in self.XY:
                row = ["%0.4f"%x,
                         "%0.4f"%self.XY[x]['hist']['skew'],
                         "%0.4f"%self.XY[x]['hist']['kurtosis'],
                         "%0.4f"%self.XY[x]['hist']['skewstat'],
                         "%0.4f"%self.XY[x]['hist']['skewpval'],
                         "%0.4f"%self.XY[x]['hist']['kurtstat'],
                         "%0.4f"%self.XY[x]['hist']['kurtpval']]
                writer.writerow(row)

        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Gmean.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = ["Voltage", "Geometric Mean", "Std Deviation"]
            writer.writerow(headers)
            for x in self.XY:
                row = ["%0.4f"%x,
                         "%0.4f"%self.XY[x]['hist']['Gmean'],
                         "%0.4f"%self.XY[x]['hist']['Gstd']]
                writer.writerow(row)

        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_RHistograms.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = []
            if self.opts.logr:
                for x in self.XY:
                    headers += ["log |R| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x]
            else:
                for x in self.XY:
                    headers += ["|R| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x]
            writer.writerow(headers)
            for i in range(0, len( self.XY[list(self.XY)[0]]['R']['hist']['bin'] ) ):
                row = []
                for x in self.XY:
                    row += ["%0.4f"%self.XY[x]['R']['hist']['bin'][i],
                    "%d"%self.XY[x]['R']['hist']['freq'][i],
                    "%0.4f"%self.XY[x]['R']['hist']['fit'][i]]
                writer.writerow(row)

        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_RHistograms_stats.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = ["Voltage", "Skew", "Kurtosis", "Skew zscore", "Skew pvalue", "Kurtosis zscore", "Kurtosis pvalue"]
            writer.writerow(headers)
            for x in self.XY:
                row = ["%0.4f"%x,
                         "%0.4f"%self.XY[x]['R']['hist']['skew'],
                         "%0.4f"%self.XY[x]['R']['hist']['kurtosis'],
                         "%0.4f"%self.XY[x]['R']['hist']['skewstat'],
                         "%0.4f"%self.XY[x]['R']['hist']['skewpval'],
                         "%0.4f"%self.XY[x]['R']['hist']['kurtstat'],
                         "%0.4f"%self.XY[x]['R']['hist']['kurtpval']]
                writer.writerow(row)


        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_LogdJdVHistograms.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = []
            for x in self.XY:
                headers += ["log |dJ/dV| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x]
            writer.writerow(headers)
            for i in range(0, len( self.GHists[list(self.GHists.keys())[0]]['hist']['bin'] ) ):
                row = []
                for x in self.GHists:
                    row += ["%0.4f"%self.GHists[x]['hist']['bin'][i],
                        "%s"%self.GHists[x]['hist']['freq'][i],
                        "%0.4f"%self.GHists[x]['hist']['fit'][i]]
                writer.writerow(row)

        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_LogdJdVHistograms_stats.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = ["Voltage", "Skew", "Kurtosis", "Skew zscore", "Skew pvalue", "Kurtosis zscore", "Kurtosis pvalue"]
            writer.writerow(headers)
            for x in self.GHists:
                row = ["%0.4f"%x,
                         "%0.4f"%self.GHists[x]['hist']['skew'],
                         "%0.4f"%self.GHists[x]['hist']['kurtosis'],
                         "%0.4f"%self.GHists[x]['hist']['skewstat'],
                         "%0.4f"%self.GHists[x]['hist']['skewpval'],
                         "%0.4f"%self.GHists[x]['hist']['kurtstat'],
                         "%0.4f"%self.GHists[x]['hist']['kurtpval']]
                writer.writerow(row)


    def WriteFilteredHistograms(self):
        '''Write the underlying histograms and associated statistics used to compute the
           Gaussian mean and variance from filtered input data according to the filtering
           cutoffs specified by the user.'''
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_filteredHistograms.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = []
            #for x in self.XY: headers += ["Log |J| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x, \
            #        "Skew (%0.4f)"%x, "Kurtosis (%0.4f)"%x, "Skew test (%0.4f)"%x, "Skew pvalue (%0.4f)"%x]
            for x in self.XY:
                headers += ["Log |J| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x]
            writer.writerow(headers)
            for i in range(0, len( self.XY[list(self.XY.keys())[0]]['filtered_hist']['bin'] ) ):
                row = []
                for x in self.XY:
                    row += ["%0.4f"%self.XY[x]['filtered_hist']['bin'][i],
                        "%s"%self.XY[x]['filtered_hist']['freq'][i],
                        "%0.4f"%self.XY[x]['filtered_hist']['fit'][i]]
                writer.writerow(row)

        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Histograms_stats.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = ["Voltage", "Skew", "Kurtosis", "Skew zscore", "Skew pvalue", "Kurtosis zscore", "Kurtosis pvalue"]
            writer.writerow(headers)
            for x in self.XY:
                row = ["%0.4f"%x,
                         "%0.4f"%self.XY[x]['filtered_hist']['skew'],
                         "%0.4f"%self.XY[x]['filtered_hist']['kurtosis'],
                         "%0.4f"%self.XY[x]['filtered_hist']['skewstat'],
                         "%0.4f"%self.XY[x]['filtered_hist']['skewpval'],
                         "%0.4f"%self.XY[x]['filtered_hist']['kurtstat'],
                         "%0.4f"%self.XY[x]['filtered_hist']['kurtpval']]
                writer.writerow(row)

    def WriteSegmentedGauss(self):
        '''Write histograms of values of J broken out by segment to catch
        hysteretic behavior without smearing it out.'''

        if not self.segments:
            logger.warning("No segments found.")
            return
        for segment in self.segments:
            rows = {}
            _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Gauss_Segment_%s.txt" % str(segment+1))
            with open(_fn, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, dialect='JV')
                headers = ["Potential (V)"]
                _maxtrace = 0
                for trace in self.segments[segment]:
                    #TODO: Fix this hack
                    if trace == 'combined':
                        continue
                    _maxtrace += 1
                    headers += ["Log|J|",
                        "Standard Deviation",
                        "Standard Error of the Mean",
                        "%s%% confidence interval" % (100*(1-self.opts.alpha)) ]
                    for x in self.segments[segment][trace]:
                        _hist = self.segments[segment][trace][x]
                        if x not in rows:
                            rows[x] = []
                        rows[x].append("%0.4f"%_hist['mean'])
                        rows[x].append("%0.4f"%_hist['std'])
                        _sem = float(_hist['std'])/np.sqrt(len(self.opts.in_files)-1 or 1)
                        rows[x].append("%0.4f"%_sem)
                        _t_val = _sem * stdtrit( len(self.opts.in_files)-1 or 1, 1 - self.opts.alpha )
                        rows[x].append("%0.4f"% _t_val)

                writer.writerow(headers)
                _V = list(rows.keys())
                _V.sort()
                for x in _V:
                    while len(rows[x]) < _maxtrace * 3:
                        rows[x] += ['-','-','-']
                        logger.warning('Filling columns for segment %i, V=%s to match %s traces.', segment, x, _maxtrace)
                    writer.writerow(["%0.4f"%x]+rows[x])

        #TODO: Don't just repeat the whole code block
        for segment in self.segments:
            rows = {}
            _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Gauss_Segments_Combined_%s.txt" % str(segment+1))
            with open(_fn, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, dialect='JV')
                headers = ["Potential (V)",
                    "Log|J|",
                    "Standard Deviation",
                    "Standard Error of the Mean",
                    "%s%% confidence interval" % (100*(1-self.opts.alpha)) ]
                for x in self.segments[segment]['combined']:
                    _hist = self.segments[segment]['combined'][x]
                    if x not in rows:
                        rows[x] = []
                    rows[x].append("%0.4f"%_hist['mean'])
                    rows[x].append("%0.4f"%_hist['std'])
                    _sem = float(_hist['std'])/np.sqrt(len(self.opts.in_files)-1 or 1)
                    rows[x].append("%0.4f"%_sem)
                    _t_val = _sem * stdtrit( len(self.opts.in_files)-1 or 1, 1 - self.opts.alpha )
                    rows[x].append("%0.4f"% _t_val)

                writer.writerow(headers)
                _V = list(rows.keys())
                _V.sort()
                for x in _V:
                    writer.writerow(["%0.4f"%x]+rows[x])

    def WriteVtrans(self):
        '''Write the Vtrans data and associated statistics.'''
        for key in ('pos', 'neg'):
            if key not in self.FN:
                logger.warning("%s not found in Fowler Nordheim data, skipping output." , key)
                continue
            _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Vtrans_"+key+".txt")
            with open(_fn, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, dialect='JV')
                writer.writerow(["Vtrans (eV)","Frequency",
                    "Gauss Fit (mean: %0.4f, Standard Deviation: %f)"%(self.FN[key]['mean'], self.FN[key]['std'])])
                data = {}
                for i in range(0, len(self.FN[key]['bin'])):
                    data[self.FN[key]['bin'][i]] = (self.FN[key]['freq'][i],self.FN[key]['fit'][i])
                for x in sorted(data.keys()):
                    writer.writerow(['%0.4f'%x,'%d'%data[x][0],'%0.2f'%data[x][1]])

        _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Vtrans_stats.txt")
        with open(_fn, 'w') as _fh:
            for key in ('pos', 'neg'):
                if key not in self.FN:
                    logger.warning("%s not found in Fowler Nordheim data, skipping output." , key)
                    continue
                _fh.write('--- %s ---\n' % key)
                _fh.write('Skew: %s\n' % self.FN[key]['skew'])
                _fh.write('Skew z-score: %s\n' % self.FN[key]['skewstat'])
                _fh.write('Skew p-val: %s\n' % self.FN[key]['skewpval'])
                _fh.write('Kurtosis: %s\n' % self.FN[key]['kurtosis'])
                _fh.write('Kurtosis z-score test: %s\n' % self.FN[key]['kurtstat'])
                _fh.write('Kurtosis p-val: %s\n' % self.FN[key]['kurtpval'])


    def WriteFN(self):
        '''Write Fowler-Nordheim plots of input data.'''
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_FN.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(['1/V'] + ['Y_%d'%x for x in range(1,len( self.XY[list(self.XY.keys())[0]]['FN'] )+1)])
            for x in self.XY:
                if x == 0.0:
                    continue
                writer.writerow(['%0.4f' % (1/x)] + list(self.XY[x]['FN']))

    def WriteFilteredGauss(self):
        '''Write the Gaussian-derived J/V data derived from the filtered input data
           according to the cutoffs specified by the user.'''
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_filteredGauss.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)",
                            "Log|J|",
                            "Standard Deviation",
                            "Standard Error of the Mean",
                            "%s%% confidence interval" % (100*(1-self.opts.alpha))])
            #Y = []
            #Yerr = []
            for x in self.XY:
                _sem = self.XY[x]['filtered_hist']['std']/np.sqrt(len(self.opts.in_files)-1 or 1)
                writer.writerow(['%f'%x,
                                '%0.4f'%self.XY[x]['filtered_hist']['mean'],
                                '%0.4f'%self.XY[x]['filtered_hist']['std'],
                                '%0.4f'% _sem,
                                '%0.4f'% (_sem * stdtrit( len(self.opts.in_files)-1 or 1, 1 - self.opts.alpha )) ])

    def WriteGauss(self):
        '''Write the Gaussian-derived data for J, R and the differential conductance data.'''
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Gauss.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)",
                            "Log|J|",
                            "Standard Deviation",
                            "Standard Error of the Mean",
                            "%s%% confidence interval" % (100*(1 - self.opts.alpha)) ])
            #Y = []
            #Yerr = []
            for x in self.XY:
                _sem = self.XY[x]['hist']['std']/np.sqrt(len(self.opts.in_files)-1 or 1)
                writer.writerow([
                        '%0.4f'%x,
                        '%0.4f'%self.XY[x]['hist']['mean'],
                        '%0.4f'%self.XY[x]['hist']['std'],
                        '%0.4f'% _sem,
                        '%0.4f'% (_sem * stdtrit( len(self.opts.in_files)-1 or 1, 1 - self.opts.alpha )) ])
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_RGauss.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            if self.opts.logr:
                writer.writerow(["Potential (V)",
                                "log |R|",
                                "Standard Deviation",
                                "%s%% confidence interval" % (100*(1 - self.opts.alpha))])
            else:
                writer.writerow(["Potential (V)",
                                "|R|",
                                "Standard Deviation",
                                "%s%% confidence interval" % (100*(1 - self.opts.alpha))])
            for x in self.XY:
                _sem = float(self.XY[x]['R']['hist']['std'])/np.sqrt(len(self.opts.in_files)-1 or 1)
                _t_val = _sem * stdtrit( len(self.opts.in_files)-1 or 1, 1 - self.opts.alpha )
                writer.writerow(['%0.4f'%x,
                                 '%0.4f'%self.XY[x]['R']['hist']['mean'],
                                 '%0.4f'%self.XY[x]['R']['hist']['std'],
                                 "%0.4f"% _t_val])
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_logdJdVGauss.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)",
                             "Log|dJ/dV|",
                            "Standard Deviation",
                            "%s%% confidence interval" % (100*(1 - self.opts.alpha))])
            for x in self.GHists:
                _sem = float(self.GHists[x]['hist']['std'])/np.sqrt(len(self.opts.in_files)-1 or 1)
                _t_val = _sem * stdtrit( len(self.opts.in_files)-1 or 1, 1 - self.opts.alpha )
                writer.writerow(['%0.4f'%x,
                                 '%0.4f'%self.GHists[x]['hist']['mean'],
                                 '%0.4f'%self.GHists[x]['hist']['std'],
                                 "%0.4f"% _t_val])
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_NDCGauss.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)",
                             "dJ/dV * V/J",
                             "Standard Deviation",
                             "%s%% confidence interval" % (100*(1 - self.opts.alpha))])
            for x in self.GHists:
                _sem = float(self.NDCHists[x]['hist']['std'])/np.sqrt(len(self.opts.in_files)-1 or 1)
                _t_val = _sem * stdtrit( len(self.opts.in_files)-1 or 1, 1 - self.opts.alpha)
                writer.writerow(['%0.4f'%x,
                                 '%0.4f'%self.NDCHists[x]['hist']['mean'],
                                 '%0.4f'%self.NDCHists[x]['hist']['std'],
                                 "%0.4f"% _t_val])

    def WriteVT(self):
        '''Write the transition voltage data plotted against potential (not 1/V).'''
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_VT.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)","|V^2/J|"])
            for x in self.XY:
                writer.writerow(['%f'%x,'%f'%self.XY[x]['VT']])

    def WriteData(self, log=False):
        '''Write LogJ or LogY (where Y is the generic Y-axis data) plotted against potential.'''
        if log:
            key,label ='LogY','LogJ'
        else:
            key, label ='Y','J'
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt")
        with open(_fn,'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)"] + ['Y_%d'%x for x in range(1,len(self.XY[list(self.XY.keys())[0]][key] )+1)])
            for x in self.XY:
                writer.writerow(["%0.4f"%x]+list(self.XY[x][key]))

    def WriteDJDV(self):
        '''Write the derivative dJ/dV plotted against potential.'''
        label = 'DJDV'
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt")
        with open(_fn,'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)"] + ['DJDV_%d'%x for x in range(1,len(self.DJDV[list(self.DJDV.keys())[0]])+1)])
            X = list(self.DJDV.keys())
            X.sort()
            for x in X:
                writer.writerow(["%0.4f"%x]+self.DJDV[x])

    def WriteNDC(self):
        '''Write the normalized differential conductance plotted against potenial.'''
        label = 'NDC'
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt")
        with open(_fn,'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)"] + ['NDC_%d'%x for x in range(1,len(self.DJDV[list(self.DJDV.keys())[0]])+1)])
            X = list(self.NDC.keys())
            X.sort()
            for x in X:
                writer.writerow(["%0.4f"%x]+self.NDC[x])

    def WriteFiltered(self):
        '''Write the filtered J/V data using the input cutoffs provided by the user.'''
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_filtered.txt")
        with open(_fn,'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            for l in self.filtered:
                writer.writerow(l)

    def WriteLag(self):
        '''Write the lag plots of the J/V data.'''
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_lag.txt")
        lenrow = len( self.XY[list(self.XY)[0]]['lag'][0] )
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            labels = []
            for x in self.XY:
                labels += ["Log|J1| @ %0.2fV"%x,"Log|J2| @ %0.2fV"%x]
            writer.writerow(labels)
            for i in range(0,lenrow):
                row = []
                for x in self.XY:
                    try:
                        row.append(self.XY[x]['lag'][0][i])
                        row.append(self.XY[x]['lag'][1][i])
                    except IndexError:
                        continue
                writer.writerow(row)

    def WriteRData(self):
        '''Write the rectification plotted against potential.'''
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Rdata.txt")
        with open(_fn,'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Potential (V)"] + ['R_%d'%x for x in range(1,len(self.XY[list(self.XY)[0]]['R'] )+1)])
            for x in self.XY:
                writer.writerow(["%0.4f"%x]+list(self.XY[x]['R']))

    def WriteGHistogram(self):
        '''Write a contour plot of conductance.'''
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_GHistogram.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = ["Potential (V)", "Log dJ/dV", "Frequency"]
            writer.writerow(headers)

            for x in self.GHists:
                for i in range(0, len(self.GHists[x]['hist']['bin'])):
                    row = [x,"%0.4f"%self.GHists[x]['hist']['bin'][i],"%s"%self.GHists[x]['hist']['freq'][i]]
                    writer.writerow(row)
                writer.writerow([])


    def WriteGMatrix(self, label):
        '''Write a matlab-style colormap maxtrix of G or NDC.'''
        if label == 'GMatrix':
            Hists = self.GHists
        elif label == 'NDCMatrix':
            Hists = self.NDCHists
        else:
            return

        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_%s.txt" % label)
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            x,y,z = [],[],[]
            for i in range(0, len(Hists[list(Hists.keys())[0]]['hist']['bin'])):
                for v in Hists:
                    x.append(v)
                    y.append(Hists[v]['hist']['bin'][i])
                    z.append(Hists[v]['hist']['freq'][i])
            x,y,z = np.array(x),np.array(y),np.array(z)
            xmin,xmax = x.min(),x.max()
            if label == 'GMatrix':
                ymin,ymax = self.opts.mlow, self.opts.mhi
            elif label == 'NDCMatrix':
                ymin,ymax = self.opts.ndc_mlow, self.opts.ndc_mhi
            else:
                ymin,ymax = y.min(),y.max()
            xi=np.linspace(xmin,xmax,200)
            yi=np.linspace(ymin,ymax,200)
            X,Y= np.meshgrid(xi,yi)
            Z = griddata((x, y), z, (X, Y),fill_value=0,method='cubic')

            for _t in enumerate(Z):
                zi = []
                for z in _t[1]:
                    if z < 0:
                        zi.append(0)
                    else:
                        zi.append(z)
                writer.writerow(zi)
        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_%s_Labels.txt" % label)
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            headers = []
            for x in X[0]:
                headers += ["%0.4f"%x]
            writer.writerow(headers)
            headers = []
            for _t in enumerate(Y):
                headers += ['%0.4f'%_t[1][0]]
            writer.writerow(headers)

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

        _fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_"+bfn+".txt")
        with open(_fn , 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            if len(labels):
                writer.writerow(labels)

            for n in range(0, lencola):
                row = []
                for _t in enumerate(dataset):
                    row.append(_t[1][n])
                writer.writerow(row)


    def WriteGNUplot(self, gpinbn, tocopy=None):
        '''Write GNUPlot sciprts for plotting the ASCII files output by the other writers.'''
        tocopy = tocopy or []
        absdir = os.path.dirname(os.path.abspath(__file__))
        tdir = os.path.join(absdir,'../../templates/')
        gpintp = os.path.join(tdir,gpinbn+'.gpin')
        try:
            nsub = str(open(gpintp,'rt').read()).count('%s')
        except FileNotFoundError:
            logger.warning("Could not read template file %s" , gpintp)
            return
        ssub = []
        for _ in range(nsub):
            ssub.append(self.opts.outfile)
        txt = open(gpintp, 'rt').read() % tuple(ssub)
        _fh = open(os.path.join(self.opts.out_dir,self.opts.outfile+'_'+gpinbn+'.gpin'), 'wt')
        _fh.write(txt)
        _fh.close()
        for _fn in tocopy:
            copyfile(os.path.join(tdir,_fn), \
                    os.path.join(self.opts.out_dir,_fn))
