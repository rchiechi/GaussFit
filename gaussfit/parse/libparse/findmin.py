from gaussfit.parse.libparse.util import throwimportwarning
from gaussfit.parse.libparse.dohistogram import dohistogram
from SLM.extract import findvtrans

try:
    import numpy as np
    import scipy.interpolate
    import pandas as pd
except ImportError as msg:
    throwimportwarning(msg)


def findmin(self):
    '''
    Find the troughs of ln(Y/X^2) vs. 1/X plots
    i.e., Vtrans, by either interpolating the data with
    a spline function and finding X where dY/dX = 0
    that gives the most negative value of Y (opts.smooth)
    or simply the most negative value of Y (! opts.smooth)
    '''
    self.logger.info("* * * * * * Computing Vtrans * * * * * * * *")
    self.loghandler.flush()
    neg_min_x, pos_min_x = [], []
    # vmin, vmax = self.df.V.min(), self.df.V.max()
    # xneg, xpos = np.linspace(vmin, 0, 250), np.linspace(vmax, 0, 250)
    tossed = 0
    if self.opts.skipohmic:
        # Vtrans has no physical meaning for curves with negative derivatives
        self.logger.info("Skipping %s (out of %s) non-tunneling traces for Vtrans calculation.",
                         len(self.ohmic), (len(self.avg.index.levels[0])))

    for trace in self.avg.index.levels[0]:
        self.SLM['Vtpos'][trace] = np.nan
        self.SLM['Vtneg'][trace] = np.nan
        if self.opts.skipohmic and trace in self.ohmic:
            tossed += 1
            continue
        FN = findvtrans(self.avg.loc[trace].index.values, self.avg.loc[trace]['J'].values, logger=self.logger)
        if not FN['err']:
            neg_min_x.append(FN['vt_neg'])
            pos_min_x.append(FN['vt_pos'])
            self.SLM['Vtpos'][trace] = FN['vt_neg']
            self.SLM['Vtneg'][trace] = FN['vt_pos']
    if tossed:
        self.logger.warning("Tossed %d compliance traces during FN calculation.", tossed)
    neg_min_x = np.array(list(filter(np.isfinite, neg_min_x)))
    pos_min_x = np.array(list(filter(np.isfinite, pos_min_x)))
    if not len(neg_min_x) or not len(pos_min_x):
        self.logger.error("Did not parse any FN values!")
        self.FN["neg"], self.FN["pos"] = {}, {}
    else:
        self.FN["neg"] = dohistogram(self.logqueue, neg_min_x, label="Vtrans(-)", density=True)
        self.FN["pos"] = dohistogram(self.logqueue, pos_min_x, label="Vtrans(+)", density=True)


def old_findmin(self):
    '''
    Find the troughs of ln(Y/X^2) vs. 1/X plots
    i.e., Vtrans, by either interpolating the data with
    a spline function and finding X where dY/dX = 0
    that gives the most negative value of Y (opts.smooth)
    or simply the most negative value of Y (! opts.smooth)
    '''
    self.logger.info("* * * * * * Computing Vtrans * * * * * * * *")
    self.loghandler.flush()
    neg_min_x, pos_min_x = [], []
    vmin, vmax = self.df.V.min(), self.df.V.max()
    xneg, xpos = np.linspace(vmin, 0, 250), np.linspace(vmax, 0, 250)
    tossed = 0
    if self.opts.skipohmic:
        # Vtrans has no physical meaning for curves with negative derivatives
        self.logger.info("Skipping %s (out of %s) non-tunneling traces for Vtrans calculation.",
                         len(self.ohmic), (len(self.avg.index.levels[0])))

    for trace in self.avg.index.levels[0]:
        self.SLM['Vtpos'][trace] = np.nan
        self.SLM['Vtneg'][trace] = np.nan
        if self.opts.skipohmic and trace in self.ohmic:
            tossed += 1
            continue
        try:
            err = [False, False]
            self.logger.debug('Finding minimum FN plot from derivative.')
            try:
                splpos = scipy.interpolate.UnivariateSpline(np.array(
                                                            self.avg.loc[trace].index[self.avg.loc[trace].index > 0]),
                                                            self.avg.loc[trace]['FN'][self.avg.loc[trace].index > 0].values,
                                                            k=4)
                pos_min_x += list(np.array(splpos.derivative().roots()))
            except Exception as msg:
                self.logger.warning('Error finding derivative of FN(+), falling back to linear interpolation. %s', str(msg))
                err[0] = True
            try:
                splneg = scipy.interpolate.UnivariateSpline(np.array(
                                                            self.avg.loc[trace].index[self.avg.loc[trace].index < 0]),
                                                            self.avg.loc[trace]['FN'][self.avg.loc[trace].index < 0].values,
                                                            k=4)
                neg_min_x += list(np.array(splneg.derivative().roots()))
            except Exception as msg:
                self.logger.warning('Error finding derivative of FN(–), falling back to linear interpolation. %s', str(msg))
                err[1] = True
            if err == (False, False):
                continue

            self.logger.debug('Finding minimum of interpolated FN plot.')
            splpos = scipy.interpolate.interp1d(np.array(self.avg.loc[trace].index[self.avg.loc[trace].index > 0]),
                                                self.avg.loc[trace]['FN'][self.avg.loc[trace].index > 0].values,
                                                kind='linear', fill_value='extrapolate')
            splneg = scipy.interpolate.interp1d(np.array(self.avg.loc[trace].index[self.avg.loc[trace].index < 0]),
                                                self.avg.loc[trace]['FN'][self.avg.loc[trace].index < 0].values,
                                                kind='linear', fill_value='extrapolate')
            xy = {'X': [], 'Y': []}
            for x in xneg:
                if not np.isfinite(x):
                    continue
                xy['Y'].append(float(splneg(x)))
                xy['X'].append(float(x))
            for x in xpos:
                if not np.isfinite(x):
                    continue
                xy['Y'].append(float(splpos(x)))
                xy['X'].append(float(x))
            fndf = pd.DataFrame(xy)
            pidx = fndf[fndf.X > 0]['Y'].idxmin()
            nidx = fndf[fndf.X < 0]['Y'].idxmin()
            if err[0]:
                try:
                    splpos = scipy.interpolate.UnivariateSpline(
                        fndf['X'][pidx-20:pidx+20].values, fndf['Y'][pidx-20:pidx+20].values, k=4)
                    pos_min_x += list(np.array(splpos.derivative().roots()))
                except Exception as msg:  # pylint: disable=broad-except
                    # TODO: Figure out how to catch all the scipy errors
                    self.logger.warning('Error finding FN(+) minimum from interpolated derivative, falling back to minimum. %s', str(msg))
                    pos_min_x.append(np.mean(fndf['X'][pidx-20:pidx+20].values))
            if err[1]:
                try:
                    splneg = scipy.interpolate.UnivariateSpline(fndf['X'][nidx-20:nidx+20].values,
                                                                fndf['Y'][nidx-20:nidx+20].values, k=4)
                    neg_min_x += list(np.array(splneg.derivative().roots()))
                except Exception as msg:  # pylint: disable=broad-except
                    # TODO: Figure out how to catch all the scipy errors
                    self.logger.warning('Error finding FN(–) minimum from interpolated derivative, falling back to minimum. %s', str(msg))
                    neg_min_x.append(np.mean(fndf['X'][nidx-20:nidx+20].values))

            self.SLM['Vtpos'][trace] = pos_min_x[-1]
            self.SLM['Vtneg'][trace] = neg_min_x[-1]
        except (IndexError, ValueError) as msg:
            self.logger.warning("Error finding minimum in trace: %s", str(msg))

    if tossed:
        self.logger.warning("Tossed %d compliance traces during FN calculation.", tossed)
    neg_min_x = np.array(list(filter(np.isfinite, neg_min_x)))
    pos_min_x = np.array(list(filter(np.isfinite, pos_min_x)))
    if not len(neg_min_x) or not len(pos_min_x):
        self.logger.error("Did not parse any FN values!")
        self.FN["neg"], self.FN["pos"] = {}, {}
    else:
        self.FN["neg"] = dohistogram(self.logqueue, neg_min_x, label="Vtrans(-)", density=True)
        self.FN["pos"] = dohistogram(self.logqueue, pos_min_x, label="Vtrans(+)", density=True)
