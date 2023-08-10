import os
import csv
import logging

logger = logging.getLogger('output')


def WriteData(self, log=False):
    '''Write LogJ or LogY (where Y is the generic Y-axis data) plotted against potential.'''
    if log:
        key, label = 'LogY', 'LogJ'
    else:
        key, label = 'Y', 'J'
    _fn = os.path.join(self.opts.out_dir, f"{self.opts.outfile}_{label}.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Potential (V)"] + ['Y_%d' % x for x in range(1, len(self.XY[list(self.XY.keys())[0]][key]) + 1)])
        for x in self.XY:
            writer.writerow(["%0.4f" % x] + list(self.XY[x][key]))

    key += '_nofirst'
    _fn = os.path.join(self.opts.out_dir, f"{self.opts.outfile}_{label}_noFirst.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Potential (V)"] + ['Y_%d' % x for x in range(1, len(self.XY[list(self.XY.keys())[0]][key]) + 1)])
        for x in self.XY:
            writer.writerow(["%0.4f" % x] + list(self.XY[x][key]))
