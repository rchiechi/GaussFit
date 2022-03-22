import os
import csv
import logging

logger = logging.getLogger('output')


def WriteFiltered(self):
    '''Write the filtered J/V data using the input cutoffs provided by the user.'''
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_filtered.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        for _l in self.filtered:
            writer.writerow(_l)
