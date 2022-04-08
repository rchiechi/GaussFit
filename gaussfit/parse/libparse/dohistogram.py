import logging
from logging.handlers import QueueHandler
from gaussfit.parse.libparse.util import throwimportwarning, signedgmean, lorenz, gauss
from gaussfit.parse.libparse import dohistogram
try:
    import numpy as np
    from scipy.optimize import curve_fit
    from scipy.stats import skew, skewtest, kurtosis, kurtosistest
except ImportError as msg:
    throwimportwarning(msg)


def dohistogram(que, opts, Y, **kwargs):
    # label="", density=False):
    '''
    Return a histogram of Y-values and a gaussian
    fit of the histogram, excluding values that
    exceed either the compliance limit (for current
    or current-density) or the ceiling for R. We
    would like to include all data in the histogram,
    but outliers sometimes confuse the fitting
    routine, which defeats the purpose of machine-fitting
    '''

    defaultKwargs = {'label': '', 'density': False}
    kwargs = {**defaultKwargs, **kwargs}
    logger = logging.getLogger(__package__+".dohistogram")
    logger.addHandler(QueueHandler(que))

    def __handlematherror(msg):
        # TODO we can now split out the file name with the bad data in it!
        logger.warning("Encountered this error while constructing histogram: %s", str(msg), exc_info=False)
        bins = np.array([0., 0., 0., 0.])
        freq = np.array([0., 0., 0., 0.])
        return bins, freq

    try:
        yrange = (Y.min(), Y.max())
    except ValueError as msg:
        logger.error("Error ranging data for histogram: %s", str(msg))
        yrange = (0, 0)

    if kwargs['label'] == "J" or kwargs['label'] == "lag":
        Y = Y[Y <= opts.compliance]
        if yrange != (0, 0):
            yrange = (Y.min()-1, Y.max()+1)
    if kwargs['label'] == "R":
        Y = Y[Y <= opts.maxr]
    if kwargs['label'] in ('DJDV', 'NDC'):
        nbins = opts.heatmapbins
    else:
        nbins = opts.bins
    if len(Y) < 10:
        logger.warning("Histogram with only %d points.", len(Y))
    try:
        # TODO Why not offer density plots as an option?
        freq, bins = np.histogram(Y, range=yrange, bins=nbins, density=kwargs['density'])
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
        if opts.lorenzian:
            coeff, covar = curve_fit(lorenz, bin_centers, freq, p0=p0, maxfev=opts.maxfev)
            hist_fit = lorenz(bin_centers, *coeff)
        else:
            coeff, covar = curve_fit(gauss, bin_centers, freq, p0=p0, maxfev=opts.maxfev)
            hist_fit = gauss(bin_centers, *coeff)
    except RuntimeError:
        if opts.maxfev > 100:
            logger.warning("|%s| Fit did not converge", kwargs['label'], exc_info=False)
    except ValueError as msg:
        logger.warning("|%s| Skipping data with ridiculous numbers in it (%s)", kwargs['label'], str(msg), exc_info=False)
        # coeff=p0
    except FloatingPointError as msg:
        logger.error("|%s| Encountered floating point error fitting Guasian: %s", kwargs['label'], str(msg), exc_info=False)

    try:
        skewstat, skewpval = skewtest(freq)
        kurtstat, kurtpval = kurtosistest(freq)
    except ValueError as msg:
        logger.error("|%s| Could not perform skewtest: %s", kwargs['label'], str(msg), exc_info=False)
        skewstat, skewpval, kurtstat, kurtpval = 0.0, 0.0, 0.0, 0.0
    return {"bin": bin_centers, "freq": freq, "mean": coeff[1], "std": coeff[2],
            "var": coeff[2], "bins": bins, "fit": hist_fit, "Gmean": Ym, "Gstd": Ys,
            "skew": skew(freq), "kurtosis": kurtosis(freq), "skewstat": skewstat, "skewpval": skewpval,
            "kurtstat": kurtstat, "kurtpval": kurtpval}
