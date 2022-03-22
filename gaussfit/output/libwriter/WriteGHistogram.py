import os
import csv
import logging

logger = logging.getLogger('output')


def WriteGHistogram(self):
    '''Write a contour plot of conductance.'''
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_GHistogram.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        headers = ["Potential (V)", "Log dJ/dV", "Frequency"]
        writer.writerow(headers)

        for x in self.GHists:
            for i in range(0, len(self.GHists[x]['hist']['bin'])):
                row = [x, "%0.4f" % self.GHists[x]['hist']['bin'][i], "%s" % self.GHists[x]['hist']['freq'][i]]
                writer.writerow(row)
            writer.writerow([])
