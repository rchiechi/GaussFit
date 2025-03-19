import os
import csv
import logging
import numpy as np
from scipy.special import stdtrit

logger = logging.getLogger('output')


def WriteGauss(self):
    '''Write the Gaussian-derived data for J, R and the differential conductance data.'''
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile + "_Gauss.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Potential (V)",
                         "Log|J|",
                         "Standard Deviation",
                         "Standard Error of the Mean",
                         "%s%% confidence interval" % (100 * (1 - self.opts.alpha))])

        for x in self.XY:
            _sem = self.XY[x]['hist']['std'] / np.sqrt(self.opts.degfree - 1 or 1)
            writer.writerow([
                            '%0.4f' % x,
                            '%0.4f' % self.XY[x]['hist']['mean'],
                            '%0.4f' % self.XY[x]['hist']['std'],
                            '%0.4f' % _sem,
                            '%0.4f' % (_sem * stdtrit(self.opts.degfree - 1 or 1, 1 - self.opts.alpha))])

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile + "_ln_Gauss.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Potential (V)",
                         "Ln|J|",
                         "Standard Deviation",
                         "Standard Error of the Mean",
                         "%s%% confidence interval" % (100 * (1 - self.opts.alpha))])

        for x in self.XY:
            _sem = self.XY[x]['ln_hist']['std'] / np.sqrt(self.opts.degfree - 1 or 1)
            writer.writerow([
                            '%0.4f' % x,
                            '%0.4f' % self.XY[x]['ln_hist']['mean'],
                            '%0.4f' % self.XY[x]['ln_hist']['std'],
                            '%0.4f' % _sem,
                            '%0.4f' % (_sem * stdtrit(self.opts.degfree - 1 or 1, 1 - self.opts.alpha))])

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile + "_Gauss_noFirstTraces.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Potential (V)",
                         "Log|J|",
                         "Standard Deviation",
                         "Standard Error of the Mean",
                         "%s%% confidence interval" % (100 * (1 - self.opts.alpha))])

        for x in self.XY:
            _sem = self.XY[x]['hist_nofirst']['std']/np.sqrt(self.opts.degfree - 1 or 1)
            writer.writerow([
                            '%0.4f' % x,
                            '%0.4f' % self.XY[x]['hist_nofirst']['mean'],
                            '%0.4f' % self.XY[x]['hist_nofirst']['std'],
                            '%0.4f' % _sem,
                            '%0.4f' % (_sem * stdtrit(self.opts.degfree - 1 or 1, 1 - self.opts.alpha))])

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_RGauss.txt")
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
            _sem = float(self.XY[x]['R']['hist']['std'])/np.sqrt(self.opts.degfree - 1 or 1)
            _t_val = _sem * stdtrit(self.opts.degfree - 1 or 1, 1 - self.opts.alpha)
            writer.writerow(['%0.4f' % x,
                             '%0.4f' % self.XY[x]['R']['hist']['mean'],
                             '%0.4f' % self.XY[x]['R']['hist']['std'],
                             "%0.4f" % _t_val])
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_logdJdVGauss.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Potential (V)",
                         "Log|dJ/dV|",
                         "Standard Deviation",
                         "%s%% confidence interval" % (100*(1 - self.opts.alpha))])
        for x in self.GHists:
            _sem = float(self.GHists[x]['hist']['std'])/np.sqrt(self.opts.degfree - 1 or 1)
            _t_val = _sem * stdtrit(self.opts.degfree - 1 or 1, 1 - self.opts.alpha)
            writer.writerow(['%0.4f' % x,
                             '%0.4f' % self.GHists[x]['hist']['mean'],
                             '%0.4f' % self.GHists[x]['hist']['std'],
                             "%0.4f" % _t_val])
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_NDCGauss.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Potential (V)",
                         "dJ/dV * V/J",
                         "Standard Deviation",
                         "%s%% confidence interval" % (100*(1 - self.opts.alpha))])
        for x in self.NDCHists:
            _sem = float(self.NDCHists[x]['hist']['std'])/np.sqrt(self.opts.degfree - 1 or 1)
            _t_val = _sem * stdtrit(self.opts.degfree - 1 or 1, 1 - self.opts.alpha)
            writer.writerow(['%0.4f' % x,
                             '%0.4f' % self.NDCHists[x]['hist']['mean'],
                             '%0.4f' % self.NDCHists[x]['hist']['std'],
                             "%0.4f" % _t_val])
