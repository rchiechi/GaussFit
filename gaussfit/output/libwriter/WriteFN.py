import os
import csv
import logging

logger = logging.getLogger('output')


def WriteFN(self):
    '''Write Fowler-Nordheim plots of input data.'''
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_FN.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(['1/V'] + ['Y_%d' % x for x in range(1, len(self.XY[list(self.XY.keys())[0]]['FN'])+1)])
        for x in self.XY:
            if x == 0.0:
                continue
            writer.writerow(['%0.4f' % (1/x)] + list(self.XY[x]['FN']))
