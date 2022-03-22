import os
import csv
import logging

logger = logging.getLogger('output')


def WriteRData(self):
    '''Write the rectification plotted against potential.'''
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Rdata.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Potential (V)"] + ['R_%d' % x for x in range(1, len(self.XY[list(self.XY)[0]]['R'])+1)])
        for x in self.XY:
            writer.writerow(["%0.4f" % x]+list(self.XY[x]['R']))
