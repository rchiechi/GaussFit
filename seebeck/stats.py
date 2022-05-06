import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.linear_model import SGDRegressor
from scipy.optimize import curve_fit


def linear_fit(opts, raw_data):
    title = os.path.basename(opts.in_files[0])
    idx = opts.col_to_parse
    y_label = raw_data[list(raw_data.keys())[0]]['labels'][idx]
    logging.info('[linear_fit] Parsing column labeled %s.', y_label)

    if opts.truetemp:
        x = np.array([[float(raw_data[f'DT{dt}']['dt'])]] for dt in opts.dTn)
    else:
        x = np.array([[float(x)] for x in opts.dTn])
    __x2 = []
    for dt in opts.dTn:
        __x2 += [[float(dt)]] * len(raw_data[f'DT{dt}']['data'][idx])
    x2 = np.array(__x2)
    lr = SGDRegressor(max_iter=100000)       # Gaussian mean-based LS

    lr.fit(x, [np.mean(raw_data[f'DT{dt}']['data'][idx]) for dt in opts.dTn])

    lr2 = SGDRegressor(max_iter=100000)      # All data-based LS
    __sum = []
    for dt in opts.dTn:
        __sum += raw_data[f'DT{dt}']['data'][idx]

    lr2.fit(x2, __sum)

    lr3 = SGDRegressor(loss='epsilon_insensitive', epsilon=0, max_iter=100000)       # All data-based LAD
    lr3.fit(x2, __sum)
    plt.errorbar(x.reshape(-1), [np.mean(raw_data[f'DT{dt}']['data'][idx]) for dt in opts.dTn],
                 yerr=[np.std(raw_data[f'DT{dt}']['data'][idx]) for dt in opts.dTn], fmt="o", color='black')
    plt.plot(x, lr.coef_[0] * x + lr.intercept_[0], label='MeanBasedLS')
    plt.plot(x2, lr2.coef_[0] * x2 + lr.intercept_[0], label='AllDataBasedLS')
    plt.plot(x2, lr3.coef_[0] * x2 + lr.intercept_[0], label='AllDataBasedLAD')
    plt.xlim(right=15)
    plt.xlabel('ΔT')
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend()
    plt.savefig(os.path.join(opts.out_dir, f'{title}_{y_label}_plot.png'))
    logging.info("successfully generated graph")       # log: print to terminal when its going well
    plt.clf()

    return lr.coef_[0], lr2.coef_[0], lr3.coef_[0]



def GHistograms(opts, raw_data):
    histograms = {}
    for dt in opts.dTn:
        histograms[f'DT{dt}'] = __Ghistrogram(raw_data[f'DT{dt}']['data'][opts.col_to_parse])
    x_label = raw_data[list(raw_data.keys())[0]]['labels'][opts.col_to_parse]
    plotGHist(opts, histograms, x_label)
    return histograms

def __Ghistrogram(array_of_data):
    G = np.array(array_of_data)
    n_bins = int(len(G)/10)
    freq, bins = np.histogram(G, bins=n_bins, density=True)
    bin_centers = (bins[:-1] + bins[1:])/2
    p0 = [1., np.mean(G), abs(np.std(G))]
    coeff = p0
    covar = None
    hist_fit = np.array([0]*len(bin_centers))
    try:
        coeff, covar = curve_fit(gauss, bin_centers, freq, p0=p0, maxfev=1000)
        hist_fit = gauss(bin_centers, *coeff)
    except RuntimeError as msg:
        logging.error('Fit did not converge: %s', msg)
    except ValueError as msg:
        logging.warning('Skipping ridiculous numbers: %s', msg)
    histogram = {"bins": bin_centers, "bin_edges": bins, "freq": freq, "mean": coeff[1], "std": coeff[2],
                 "var": coeff[2], "fit": hist_fit, "Gmean": np.mean(G), "Gstd": np.std(G), "covar": covar}
    return histogram


def plotGHist(opts, hist, x_label):
    title = os.path.basename(opts.in_files[0])
    colors = ['blue', 'red', 'green', 'organge', 'purple']
    i = 0
    for dt in opts.dTn:
        plt.plot(hist[f'DT{dt}']['bins'], hist[f'DT{dt}']['fit'], color=colors[i], label=f'ΔT{dt}')
        plt.hist(hist[f'DT{dt}']['bins'], weights=hist[f'DT{dt}']['freq'], bins=len(hist[f'DT{dt}']['bins']))
        if i == len(colors):
            i = 0
        else:
            i += 1
    plt.legend()
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel('counts')
    plt.savefig(os.path.join(opts.out_dir, f'{title}_{x_label}_GHistograms.png'))
    plt.clf()
    logging.info("successfully generated GHistograms")


def gauss(x, *p):
    """
    Return a gaussian function
    """
    a, mu, sigma = p
    return a*np.exp(-(x-mu)**2/(2.*sigma**2))


def plot_histograms(input_dir, output_dir, raw_data):
    title = input_dir.split("/")[-1]
    if len(raw_data) == 6:
        title += '_dT'
    for i in raw_data[:3]:
        i.sort()
    # TODO maybe just have this use gauss instead of scipy.stats.norm
    DT4pdf = norm.pdf(raw_data[0], np.mean(raw_data[0]), np.std(raw_data[0]))
    DT8pdf = norm.pdf(raw_data[1], np.mean(raw_data[1]), np.std(raw_data[1]))
    DT12pdf = norm.pdf(raw_data[2], np.mean(raw_data[2]), np.std(raw_data[2]))
    nbins = 160
    plt.hist(raw_data[0], bins=nbins, label='DT4')
    plt.plot(raw_data[0], nbins * DT4pdf, color='blue')     # the pdf and hist dont line up accurately
    plt.hist(raw_data[1], bins=nbins, label='DT8')      # the nbins multiplier is arbitrary to try and fix (doesnt work)
    plt.plot(raw_data[1], nbins * DT8pdf, color='red')
    plt.hist(raw_data[2], bins=nbins, label='DT12')
    plt.plot(raw_data[2], nbins * DT12pdf, color='green')
    # plt.xlim(left=-275, right=25)     # hist scale used by korean paper. for comparison
    plt.xlabel('dV')
    plt.ylabel('counts')
    plt.title(title)
    plt.legend()
    plt.savefig(f'{output_dir}/_{title}_dV_Histograms.png')
    plt.clf()
    print("successfully generated histograms")
