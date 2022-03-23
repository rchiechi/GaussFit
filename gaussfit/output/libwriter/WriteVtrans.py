import os
import csv
import logging

logger = logging.getLogger('output')


def WriteVtrans(self):
    '''Write the Vtrans data and associated statistics.'''
    for key in ('pos', 'neg'):
        if key not in self.FN:
            logger.warning("%s not found in Fowler Nordheim data, skipping output.", key)
            continue
        _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Vtrans_"+key+".txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow(["Vtrans (eV)", "Frequency",
                             "Gauss Fit (mean: %0.4f, Standard Deviation: %f)" % (
                                    self.FN[key]['mean'], self.FN[key]['std'])])
            data = {}
            for i in range(0, len(self.FN[key]['bin'])):
                data[self.FN[key]['bin'][i]] = (self.FN[key]['freq'][i], self.FN[key]['fit'][i])
            for x in sorted(data.keys()):
                writer.writerow(['%0.4f' % x, '%d' % data[x][0], '%0.2f' % data[x][1]])

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Vtrans_stats.txt")
    with open(_fn, 'w') as _fh:
        for key in ('pos', 'neg'):
            if key not in self.FN:
                logger.warning("%s not found in Fowler Nordheim data, skipping output.", key)
                continue
            _fh.write('--- %s ---\n' % key)
            _fh.write('Skew: %s\n' % self.FN[key]['skew'])
            _fh.write('Skew z-score: %s\n' % self.FN[key]['skewstat'])
            _fh.write('Skew p-val: %s\n' % self.FN[key]['skewpval'])
            _fh.write('Kurtosis: %s\n' % self.FN[key]['kurtosis'])
            _fh.write('Kurtosis z-score test: %s\n' % self.FN[key]['kurtstat'])
            _fh.write('Kurtosis p-val: %s\n' % self.FN[key]['kurtpval'])