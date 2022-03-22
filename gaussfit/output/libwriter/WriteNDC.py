import os
import csv
import logging

logger = logging.getLogger('output')


def WriteNDC(self):
    '''Write the normalized differential conductance plotted against potenial.'''
    label = 'NDC'
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_"+label+".txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Potential (V)"] + ['NDC_%d' % x for x in range(1, len(self.DJDV[list(self.DJDV.keys())[0]])+1)])
        X = list(self.NDC.keys())
        X.sort()
        for x in X:
            writer.writerow(["%0.4f" % x]+self.NDC[x])
