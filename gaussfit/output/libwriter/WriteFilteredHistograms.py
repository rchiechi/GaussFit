import os
import csv
import logging

logger = logging.getLogger('output')


def WriteFilteredHistograms(self):
    '''Write the underlying histograms and associated statistics used to compute the
        Gaussian mean and variance from filtered input data according to the filtering
        cutoffs specified by the user.'''
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_filteredHistograms.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        headers = []
        # for x in self.XY: headers += ["Log |J| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x, \
        #        "Skew (%0.4f)"%x, "Kurtosis (%0.4f)"%x, "Skew test (%0.4f)"%x, "Skew pvalue (%0.4f)"%x]
        for x in self.XY:
            headers += ["Log |J| (%0.4f)" % x, "Frequency (%0.4f)" % x, "Fit (%0.4f)" % x]
        writer.writerow(headers)
        for i in range(0, len(self.XY[list(self.XY.keys())[0]]['filtered_hist']['bin'])):
            row = []
            for x in self.XY:
                row += ["%0.4f" % self.XY[x]['filtered_hist']['bin'][i],
                        "%s" % self.XY[x]['filtered_hist']['freq'][i],
                        "%0.4f" % self.XY[x]['filtered_hist']['fit'][i]]
            writer.writerow(row)

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Histograms_stats.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        headers = ["Voltage", "Skew", "Kurtosis", "Skew zscore", "Skew pvalue", "Kurtosis zscore", "Kurtosis pvalue"]
        writer.writerow(headers)
        for x in self.XY:
            row = ["%0.4f" % x,
                   "%0.4f" % self.XY[x]['filtered_hist']['skew'],
                   "%0.4f" % self.XY[x]['filtered_hist']['kurtosis'],
                   "%0.4f" % self.XY[x]['filtered_hist']['skewstat'],
                   "%0.4f" % self.XY[x]['filtered_hist']['skewpval'],
                   "%0.4f" % self.XY[x]['filtered_hist']['kurtstat'],
                   "%0.4f" % self.XY[x]['filtered_hist']['kurtpval']]
            writer.writerow(row)
