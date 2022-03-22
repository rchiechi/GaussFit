import os
import csv
import logging

logger = logging.getLogger('output')


def WriteLag(self):
    '''Write the lag plots of the J/V data.'''
    if self.opts.nolag:
        return
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_lag.txt")
    lenrow = len(self.XY[list(self.XY)[0]]['lag'][0])
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        labels = []
        for x in self.XY:
            labels += ["Log|J1| @ %0.2fV" % x, "Log|J2| @ %0.2fV" % x]
        writer.writerow(labels)
        for i in range(0, lenrow):
            row = []
            for x in self.XY:
                try:
                    row.append(self.XY[x]['lag'][0][i])
                    row.append(self.XY[x]['lag'][1][i])
                except IndexError:
                    continue
            writer.writerow(row)
