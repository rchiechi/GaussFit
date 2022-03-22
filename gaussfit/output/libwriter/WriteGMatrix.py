import os
import csv
import logging
import numpy as np
from scipy.interpolate import griddata

logger = logging.getLogger('output')


def WriteGMatrix(self, label):
    '''Write a matlab-style colormap maxtrix of G or NDC.'''
    if label == 'GMatrix':
        Hists = self.GHists
    elif label == 'NDCMatrix':
        Hists = self.NDCHists
    elif label == 'LogJMatrix':
        Hists = self.XY
    else:
        return

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_%s.txt" % label)
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        x, y, z = [], [], []
        for i in range(0, len(Hists[list(Hists.keys())[0]]['hist']['bin'])):
            for v in Hists:
                x.append(v)
                y.append(Hists[v]['hist']['bin'][i])
                z.append(Hists[v]['hist']['freq'][i])
        x, y, z = np.array(x), np.array(y), np.array(z)
        xmin, xmax = x.min(), x.max()
        if label == 'GMatrix':
            ymin, ymax = self.opts.mlow, self.opts.mhi
        elif label == 'NDCMatrix':
            ymin, ymax = self.opts.ndc_mlow, self.opts.ndc_mhi
        else:
            ymin, ymax = y.min(), y.max()
        xi = np.linspace(xmin, xmax, 200)
        yi = np.linspace(ymin, ymax, 200)
        X, Y = np.meshgrid(xi, yi)
        Z = griddata((x, y), z, (X, Y), fill_value=0, method='cubic')

        for _t in enumerate(Z):
            zi = []
            for z in _t[1]:
                if z < 0:
                    zi.append(0)
                else:
                    zi.append(z)
            writer.writerow(zi)
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_%s_Labels.txt" % label)
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        headers = []
        for x in X[0]:
            headers += ["%0.4f" % x]
        writer.writerow(headers)
        headers = []
        for _t in enumerate(Y):
            headers += ['%0.4f' % _t[1][0]]
        writer.writerow(headers)
