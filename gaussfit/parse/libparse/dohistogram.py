from gaussfit.parse.libparse.util import throwimportwarning, signedgmean, lorenz, gauss
try:
    import numpy as np
    from scipy.optimize import curve_fit
    from scipy.stats import skew, skewtest, kurtosis, kurtosistest
except ImportError as msg:
    throwimportwarning(msg)


def dohistogram(self, Y, label="", density=False):
    '''
    Return a histogram of Y-values and a gaussian
    fit of the histogram, excluding values that
    exceed either the compliance limit (for current
    or current-density) or the ceiling for R. We
    would like to include all data in the histogram,
    but outliers sometimes confuse the fitting
    routine, which defeats the purpose of machine-fitting
    '''

    def __handlematherror(msg):
        # TODO we can now split out the file name with the bad data in it!
        self.logger.warning("Encountered this error while constructing histogram: %s", str(msg), exc_info=False)
        bins = np.array([0., 0., 0., 0.])
        freq = np.array([0., 0., 0., 0.])
        return bins, freq

    try:
        yrange = (Y.min(), Y.max())
    except ValueError as msg:
        self.logger.error("Error ranging data for histogram: %s", str(msg))
        yrange = (0, 0)

    if label == "J" or label == "lag":
        Y = Y[Y <= self.opts.compliance]
        if yrange != (0, 0):
            yrange = (Y.min()-1, Y.max()+1)
    if label == "R":
        Y = Y[Y <= self.opts.maxr]
    if label in ('DJDV', 'NDC'):
        nbins = self.opts.heatmapbins
    else:
        nbins = self.opts.bins
    if len(Y) < 10:
        self.logger.warning("Histogram with only %d points.", len(Y))
    try:
        # TODO Why not offer density plots as an option?
        freq, bins = np.histogram(Y, range=yrange, bins=nbins, density=density)
    except ValueError as msg:
        bins, freq = __handlematherror(msg)
    except FloatingPointError as msg:
        bins, freq = __handlematherror(msg)

    if len(Y):
        Ym = signedgmean(Y)
        Ys = abs(Y.std())
    else:
        Ym, Ys = 0.0, 0.0

    p0 = [1., Ym, Ys]
    bin_centers = (bins[:-1] + bins[1:])/2
    coeff = p0
    covar = None  # pylint: disable=unused-variable
    hist_fit = np.array([x*0 for x in range(0, len(bin_centers))])
    try:
        # with self.lock:
        if self.opts.lorenzian:
            coeff, covar = curve_fit(lorenz, bin_centers, freq, p0=p0, maxfev=self.opts.maxfev)
            hist_fit = lorenz(bin_centers, *coeff)
        else:
            coeff, covar = curve_fit(gauss, bin_centers, freq, p0=p0, maxfev=self.opts.maxfev)
            hist_fit = gauss(bin_centers, *coeff)
    except RuntimeError:
        if self.opts.maxfev > 100:
            self.logger.warning("|%s| Fit did not converge", label, exc_info=False)
    except ValueError as msg:
        self.logger.warning("|%s| Skipping data with ridiculous numbers in it (%s)", label, str(msg), exc_info=False)
        # coeff=p0
    except FloatingPointError as msg:
        self.logger.error("|%s| Encountered floating point error fitting Guasian: %s", label, str(msg), exc_info=False)

    try:
        skewstat, skewpval = skewtest(freq)
        kurtstat, kurtpval = kurtosistest(freq)
    except ValueError as msg:
        self.logger.error("|%s| Could not perform skewtest: %s", label, str(msg), exc_info=False)
        skewstat, skewpval, kurtstat, kurtpval = 0.0, 0.0, 0.0, 0.0
    return {"bin": bin_centers, "freq": freq, "mean": coeff[1], "std": coeff[2],
            "var": coeff[2], "bins": bins, "fit": hist_fit, "Gmean": Ym, "Gstd": Ys,
            "skew": skew(freq), "kurtosis": kurtosis(freq), "skewstat": skewstat, "skewpval": skewpval,
            "kurtstat": kurtstat, "kurtpval": kurtpval}
