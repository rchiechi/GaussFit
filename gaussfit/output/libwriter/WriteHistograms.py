import os
import csv
import logging

logger = logging.getLogger('output')


def WriteHistograms(self):
    '''Write all of the underlying histograms and associated statistics used to compute
        Gaussian mean and variance of J, R and Vtrans from the raw input data.'''
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Histograms.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        headers = []
        for x in self.XY:
            headers += ["Log |J| (%0.4f)" % x, "Frequency (%0.4f)" % x, "Fit (%0.4f)" % x]
        writer.writerow(headers)
        for i in range(0, len(self.XY[list(self.XY.keys())[0]]['hist']['bin'])):
            row = []
            for x in self.XY:
                row += ["%0.4f" % self.XY[x]['hist']['bin'][i],
                        "%s" % self.XY[x]['hist']['freq'][i],
                        "%0.4f" % self.XY[x]['hist']['fit'][i]]
            writer.writerow(row)

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Histograms_noFirstTraces.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        headers = []
        for x in self.XY:
            headers += ["Log |J| (%0.4f)" % x, "Frequency (%0.4f)" % x, "Fit (%0.4f)" % x]
        writer.writerow(headers)
        for i in range(0, len(self.XY[list(self.XY.keys())[0]]['hist_nofirst']['bin'])):
            row = []
            for x in self.XY:
                row += ["%0.4f" % self.XY[x]['hist_nofirst']['bin'][i],
                        "%s" % self.XY[x]['hist_nofirst']['freq'][i],
                        "%0.4f" % self.XY[x]['hist_nofirst']['fit'][i]]
            writer.writerow(row)

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Histograms_stats.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        headers = ["Voltage", "Skew", "Kurtosis", "Skew zscore", "Skew pvalue", "Kurtosis zscore", "Kurtosis pvalue"]
        writer.writerow(headers)
        for x in self.XY:
            row = ["%0.4f" % x,
                   "%0.4f" % self.XY[x]['hist']['skew'],
                   "%0.4f" % self.XY[x]['hist']['kurtosis'],
                   "%0.4f" % self.XY[x]['hist']['skewstat'],
                   "%0.4f" % self.XY[x]['hist']['skewpval'],
                   "%0.4f" % self.XY[x]['hist']['kurtstat'],
                   "%0.4f" % self.XY[x]['hist']['kurtpval']]
            writer.writerow(row)

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Gmean.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        headers = ["Voltage", "Geometric Mean", "Std Deviation"]
        writer.writerow(headers)
        for x in self.XY:
            row = ["%0.4f" % x,
                   "%0.4f" % self.XY[x]['hist']['Gmean'],
                   "%0.4f" % self.XY[x]['hist']['Gstd']]
            writer.writerow(row)

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_RHistograms.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        headers = []
        if self.opts.logr:
            for x in self.XY:
                headers += ["log |R| (%0.4f)" % x, "Frequency (%0.4f)" % x, "Fit (%0.4f)" % x]
        else:
            for x in self.XY:
                headers += ["|R| (%0.4f)" % x, "Frequency (%0.4f)" % x, "Fit (%0.4f)" % x]
        writer.writerow(headers)
        for i in range(0, len(self.XY[list(self.XY)[0]]['R']['hist']['bin'])):
            row = []
            for x in self.XY:
                row += ["%0.4f" % self.XY[x]['R']['hist']['bin'][i],
                        "%d" % self.XY[x]['R']['hist']['freq'][i],
                        "%0.4f" % self.XY[x]['R']['hist']['fit'][i]]
            writer.writerow(row)

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_RHistograms_stats.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        headers = ["Voltage", "Skew", "Kurtosis", "Skew zscore", "Skew pvalue", "Kurtosis zscore", "Kurtosis pvalue"]
        writer.writerow(headers)
        for x in self.XY:
            row = ["%0.4f" % x,
                   "%0.4f" % self.XY[x]['R']['hist']['skew'],
                   "%0.4f" % self.XY[x]['R']['hist']['kurtosis'],
                   "%0.4f" % self.XY[x]['R']['hist']['skewstat'],
                   "%0.4f" % self.XY[x]['R']['hist']['skewpval'],
                   "%0.4f" % self.XY[x]['R']['hist']['kurtstat'],
                   "%0.4f" % self.XY[x]['R']['hist']['kurtpval']]
            writer.writerow(row)

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_LogdJdVHistograms.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        headers = []
        for x in self.XY:
            headers += ["log |dJ/dV| (%0.4f)" % x, "Frequency (%0.4f)" % x, "Fit (%0.4f)" % x]
        writer.writerow(headers)
        for i in range(0, len(self.GHists[list(self.GHists.keys())[0]]['hist']['bin'])):
            row = []
            for x in self.GHists:
                row += ["%0.4f" % self.GHists[x]['hist']['bin'][i],
                        "%s" % self.GHists[x]['hist']['freq'][i],
                        "%0.4f" % self.GHists[x]['hist']['fit'][i]]
            writer.writerow(row)

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_LogdJdVHistograms_stats.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        headers = ["Voltage", "Skew", "Kurtosis", "Skew zscore", "Skew pvalue", "Kurtosis zscore", "Kurtosis pvalue"]
        writer.writerow(headers)
        for x in self.GHists:
            row = ["%0.4f" % x,
                   "%0.4f" % self.GHists[x]['hist']['skew'],
                   "%0.4f" % self.GHists[x]['hist']['kurtosis'],
                   "%0.4f" % self.GHists[x]['hist']['skewstat'],
                   "%0.4f" % self.GHists[x]['hist']['skewpval'],
                   "%0.4f" % self.GHists[x]['hist']['kurtstat'],
                   "%0.4f" % self.GHists[x]['hist']['kurtpval']]
            writer.writerow(row)
