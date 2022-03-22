import os
import csv
import logging

logger = logging.getLogger('output')


def WriteGeneric(self, dataset, bfn, labels=None):
    '''Write a generic set of data expecting an n-dimensional array'''
    labels = labels or []
    if labels and len(labels) != len(dataset):
        logger.error("Length of column labels does not match number of data columns for WriteGeneric!")
        return

    lencola = len(dataset[0])
    for d in dataset:
        if len(d) != lencola:
            logger.error("Length of columns differs for WriteGeneric!")
            return

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_"+bfn+".txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        if len(labels):
            writer.writerow(labels)

        for n in range(0, lencola):
            row = []
            for _t in enumerate(dataset):
                row.append(_t[1][n])
            writer.writerow(row)
