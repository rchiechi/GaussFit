import os
import csv
import logging
import numpy as np
from scipy.special import stdtrit

logger = logging.getLogger('output')


def WriteFilteredGauss(self):
    '''Write the Gaussian-derived J/V data derived from the filtered input data
        according to the cutoffs specified by the user.'''
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_filteredGauss.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Potential (V)",
                         "Log|J|",
                         "Standard Deviation",
                         "Standard Error of the Mean",
                         "%s%% confidence interval" % (100*(1-self.opts.alpha))])

        for x in self.XY:
            _sem = self.XY[x]['filtered_hist']['std']/np.sqrt(self.opts.degfree - 1 or 1)
            writer.writerow(['%f' % x,
                             '%0.4f' % self.XY[x]['filtered_hist']['mean'],
                             '%0.4f' % self.XY[x]['filtered_hist']['std'],
                             '%0.4f' % _sem,
                             '%0.4f' % (_sem * stdtrit(self.opts.degfree - 1 or 1, 1 - self.opts.alpha))])
