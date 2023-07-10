from gaussfit.parse.libparse.util import throwimportwarning

try:
    from scipy.stats import linregress
    import numpy as np
except ImportError as msg:
    throwimportwarning(msg)

def doconductance(self):
    '''
    Find the conductance using a linear regression on the first four data points.
    '''
    self.logger.info("* * * * * * Computing G * * * * * * * *")
    self.loghandler.flush()
    tossed = 0
    voltages = []
    # Make a list of V = 0 the two values of V above and below
    for _i in range(-2,3):
        _idx = self.avg.loc[0].index.tolist().index(0) + _i
        voltages.append(self.avg.loc[0].index.tolist()[_idx])
    for trace in self.avg.index.levels[0]:
        if self.opts.skipohmic and trace in self.ohmic:
            tossed += 1
            continue
        _Y = []
        for _v in voltages:
            _Y.append(self.avg.loc[trace]['J'][_v])
        _fit = linregress(voltages, _Y)
        self.logger.debug(f"G:{_fit.slope} (R={_fit.rvalue:.2f})")
        self.G[trace] = _fit.slope

    self.logger.info("Average conductance: %s", np.average(list(self.G.values())))
