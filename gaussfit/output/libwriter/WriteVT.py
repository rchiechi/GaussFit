import os
import csv
import logging

logger = logging.getLogger('output')


def WriteVT(self):
    '''Write the transition voltage data plotted against potential (not 1/V).'''
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_VT.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Potential (V)", "|V^2/J|"])
        for x in self.XY:
            writer.writerow(['%f' % x, '%f' % self.XY[x]['VT']])
