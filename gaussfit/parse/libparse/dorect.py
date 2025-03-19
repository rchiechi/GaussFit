import asyncio
from gaussfit.parse.libparse.util import throwimportwarning
from gaussfit.parse.libparse.dohistogram import dohistogram
from collections import OrderedDict

try:
    import numpy as np
except ImportError as msg:
    throwimportwarning(msg)


def dorect(self, xy):
    '''
    Divide each value of Y at +V by Y at -V
    and build a histogram of rectification, R
    also construct the unique Voltage list
    '''
    self.logger.info("* * * * * * Computing |R|  * * * * * * * * *")
    self.loghandler.flush()
    r = OrderedDict()
    R = OrderedDict()
    for x, _ in xy:
        r[x] = []
    clipped = 0
    for trace in self.avg.index.levels[0]:
        for x in self.avg.loc[trace].index[self.avg.loc[trace].index >= 0]:
            if -1*x not in self.avg.loc[trace]['J']:
                self.logger.warning("Rectification data missing voltages.")
                if self.opts.logr:
                    r[x].append(0.)
                else:
                    r[x].append(1.)
                continue
            elif x == 0.0:
                if self.opts.logr:
                    r[x].append(0.)
                else:
                    r[x].append(1.)
                continue
            if self.opts.logr:
                r[x].append(np.log10(abs(self.avg.loc[trace]['J'][x]/self.avg.loc[trace]['J'][-1*x])))
            else:
                r[x].append(abs(self.avg.loc[trace]['J'][x]/self.avg.loc[trace]['J'][-1*x]))
            if r[x][-1] > self.opts.maxr:
                clipped += 1
    for x in reversed(list(r)):
        if x >= 0:
            R[x] = {'r': np.array(r[x]), 'hist': dohistogram(np.array(r[x]), label="R", que=self.logqueue)}
            R[-1*x] = R[x]
        if x not in R:
            self.logger.warning("Unequal +/- voltages in R-plot will be filled with R=1.")
            if self.opts.logr:
                y = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])
            else:
                y = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
            R[x] = {'r': y, 'hist': dohistogram(y, label="R", que=self.logqueue)}
    if clipped:
        if self.opts.logr:
            rstr = 'log|R|'
        else:
            rstr = '|R|'
        self.logger.info("%s values of %s exceed maxR (%s)", clipped, rstr, self.opts.maxr)
    self.logger.info("R complete.")
    return R
