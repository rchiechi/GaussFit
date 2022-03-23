from gaussfit.parse.libparse.util import throwimportwarning

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
    neg_min_x, pos_min_x = [], []
    vmin, vmax = self.df.V.min(), self.df.V.max()
    xneg, xpos = np.linspace(vmin, 0, 250), np.linspace(vmax, 0, 250)
    tossed = 0
    if self.opts.skipohmic:
        # Vtrans has no physical meaning for curves with negative derivatives
        self.logger.info("Skipping %s (out of %s) non-tunneling traces for Vtrans calculation.",
                         len(self.ohmic), (len(self.avg.index.levels[0])))

    for trace in self.avg.index.levels[0]:
        if self.opts.skipohmic and trace in self.ohmic:
            tossed += 1
            continue
        try:
            if not self.opts.interpolateminfn:

                self.logger.debug('Finding FN min of plot.')
                neg_min_x.append(self.avg.loc[trace]['FN'][self.avg.loc[trace].index < 0].idxmin())
                pos_min_x.append(self.avg.loc[trace]['FN'][self.avg.loc[trace].index > 0].idxmin())

            else:
                err = [False, False]
                self.logger.debug('Finding minimum FN plot from derivative.')
                splpos = scipy.interpolate.UnivariateSpline(np.array(self.avg.loc[trace].index[self.avg.loc[trace].index > 0]),
                                                            self.avg.loc[trace]['FN'][self.avg.loc[trace].index > 0].values, k=4)
                splneg = scipy.interpolate.UnivariateSpline(np.array(self.avg.loc[trace].index[self.avg.loc[trace].index < 0]),
                                                            self.avg.loc[trace]['FN'][self.avg.loc[trace].index < 0].values, k=4)
                try:
                    pos_min_x += list(np.array(splpos.derivative().roots()))
                except ValueError as msg:
                    self.logger.warning('Error finding derivative of FN(+), falling back to linear interpolation. %s', str(msg))
                    err[0] = True
                try:
                    neg_min_x += list(np.array(splneg.derivative().roots()))
                except ValueError as msg:
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
        except ValueError as msg:
            self.logger.warning("Error finding minimum in trace: %s", str(msg))
    if tossed:
        self.logger.warning("Tossed %d compliance traces during FN calculation.", tossed)
    neg_min_x = np.array(list(filter(np.isfinite, neg_min_x)))
    pos_min_x = np.array(list(filter(np.isfinite, pos_min_x)))
    self.FN["neg"], self.FN["pos"] = self.dohistogram(neg_min_x, "Vtrans(-)", True), self.dohistogram(pos_min_x, "Vtrans(+)", True)