import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.linear_model import SGDRegressor
from scipy.optimize import curve_fit


def linear_fit(opts, raw_data, idx=1):
    title = os.path.basename(opts.in_files[0])[:-4]
    logging.info('Parsing column labeled %s.' % raw_data[list(raw_data.keys())[0]]['labels'][idx])
    if opts.truetemp:
        x = np.array([[float(raw_data[dt]['dt'])]] for dt in opts.dTn)
    else:
        x = np.array([[float(x)] for x in opts.dTn])
    __x2 = []
    for dt in opts.dTn:
        __x2 += [[float(dt)]] * len(raw_data[f'DT{dt}']['data'][idx])
    x2 = np.array(__x2)
    # TODO maybe replace SGDRegressor with own function ?
    lr = SGDRegressor(max_iter=100000)  # Gaussian mean-based LS
    lr.fit(x, [np.mean(raw_data[f'DT{dt}']['data'][idx]) for dt in opts.dTn])
    lr2 = SGDRegressor(max_iter=100000)  # All data-based LS
    __sum = []
    for dt in opts.dTn:
        __sum += raw_data[f'DT{dt}']['data'][idx]
    lr2.fit(x2, __sum)
    lr3 = SGDRegressor(loss='epsilon_insensitive', epsilon=0, max_iter=100000)  # All data-based LAD
    lr3.fit(x2, __sum)
    # TODO split this into linear_fit and Image_generator
    # generate the image:
    # plt.errorbar(x.reshape(-1), [np.mean(raw_data[f'DT{dt}']['data'][idx]) for dt in opts.dTn],
    #              yerr=[np.std(raw_data[f'DT{dt}']['data'][idx]) for dt in opts.dTn], fmt="o", color='black')
    plt.scatter(x.reshape(-1), [np.mean(raw_data[f'DT{dt}']['data'][idx]) for dt in opts.dTn],
                label='LinearFit')
    # for dt in opts.dTn:
    #     plt.scatter([float(dt)] * len(raw_data[f'DT{dt}']['data'][idx]), raw_data[f'DT{dt}']['data'][idx])
        # plt.scatter(x2[len(raw_data[0]):-len(raw_data[2])], raw_data[1])
        # plt.scatter(x2[-(len(raw_data[2])):], raw_data[2])
    plt.plot(x, lr.coef_[0] * x + lr.intercept_[0], label='MeanBasedLS')
    plt.plot(x2, lr2.coef_[0] * x2 + lr.intercept_[0], label='AllDataBasedLS')
    plt.plot(x2, lr3.coef_[0] * x2 + lr.intercept_[0], label='AllDataBasedLAD')
    plt.xlim(right=15)
    plt.xlabel('dT')
    plt.ylabel('dV')
    plt.title(title)
    plt.legend()
    plt.savefig(os.path.join(opts.out_dir, f'{title}_dVdT_plot.png'))
    logging.info("successfully generated graph")       # log: print to terminal when its going well
    plt.clf()

    return [lr.coef_[0], lr2.coef_[0], lr3.coef_[0]]


def GHistograms(*raw_data):
    G = np.array(raw_data[0])
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
        logging.warn(f'Fit did not converge: {msg}')
    except ValueError as msg:
        logging.warn(f'Skipping ridiculous numbers: {msg}')
    histogram = {"bins": bin_centers, "bin_edges": bins, "freq": freq, "mean": coeff[1], "std": coeff[2],
                 "var": coeff[2], "fit": hist_fit, "Gmean": np.mean(G), "Gstd": np.std(G)}
    return histogram


def plotGHist(input_dir, output_dir, hist):
    title = input_dir.split("/")[-1]
    # assuming hist is a list of three dicts
    plt.plot(hist[0]['bins'], hist[0]['fit'], color='blue', label='dT4')
    plt.hist(hist[0]['bins'], weights=hist[0]['freq'], bins=len(hist[0]['bins']))
    plt.plot(hist[1]['bins'], hist[1]['fit'], color='red', label='dT8')
    plt.hist(hist[1]['bins'], weights=hist[1]['freq'], bins=len(hist[1]['bins']))
    plt.plot(hist[2]['bins'], hist[2]['fit'], color='green', label='dT12')
    plt.hist(hist[2]['bins'], weights=hist[2]['freq'], bins=len(hist[2]['bins']))
    plt.legend()
    plt.title(title)
    plt.xlabel('dV')
    plt.ylabel('counts')
    plt.savefig(f'{output_dir}/_{title}_dV_GHistograms.png')
    plt.clf()
    print("successfully generated GHistograms")


def gauss(x, *p):
    """
    Return a gaussian function
    """
    a, mu, sigma = p
    return a*np.exp(-(x-mu)**2/(2.*sigma**2))


def histograms(input_dir, output_dir, raw_data):
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
