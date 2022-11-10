import os
import csv

csv.register_dialect('JV', delimiter='\t', quoting=csv.QUOTE_MINIMAL)

def doOutput(opts, **kwargs):
    base_name = os.path.join(opts.out_dir, opts.out_file)
    if 'histograms' in kwargs:
        # for dt in kwargs['histograms']:
        doHistograms(f'{base_name}_deltaVHistograms.txt', kwargs['histograms'])
    if 'linear_fit' in kwargs:
        doLinearfit(base_name, kwargs['linear_fit'])


def doHistograms(fn, histograms):
    #    histogram = {"bins": bin_centers, "bin_edges": bins, "freq": freq, "mean": coeff[1], "std": coeff[2],
    # "var": coeff[2], "fit": hist_fit, "Gmean": np.mean(G), "Gstd": np.std(G), "covar": covar}
    with open(fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        header = []
        for dt in histograms:
            header += [f'{histograms["x_label"]} (ΔT{dt})', f'Frequency (ΔT{dt})', f'Fit (ΔT{dt})']
        writer.writerow(header)
        row = []
        for dt in histograms:
            if dt == 'x_label':
                continue
            for i, _bin in enumerate(histograms[dt]['bins']):
                _freq, _fit = histograms[dt]['freq'][i], histograms[dt]['fit'][i]
                row +=         [f'{_bin:.4f}', f'{_freq:.4f}', f'{_fit:.4f}']
        writer.writerow(row)


def doLinearfit(bn, linear_fit):
    # plt.errorbar(x.reshape(-1), [np.mean(raw_data[f'DT{dt}']['data'][idx]) for dt in opts.dTn],
    #              yerr=[np.std(raw_data[f'DT{dt}']['data'][idx]) for dt in opts.dTn], fmt="o", color='black')
    # plt.plot(x, lr.coef_[0] * x + lr.intercept_[0], label='MeanBasedLS')
    # plt.plot(x2, lr2.coef_[0] * x2 + lr.intercept_[0], label='AllDataBasedLS')
    # plt.plot(x2, lr3.coef_[0] * x2 + lr.intercept_[0], label='AllDataBasedLAD')
    # return {'MeanBasedLS': (x, lr.coef_[0]),
    # 'AllDataBasedLS': (x2, lr2.coef_[0]),
    # 'AllDataBasedLAD': (x2, lr3.coef_[0]),
    # 'intercept': lr.intercept_[0],
    # 'XY': (X, Y, Y_err)}
    # fn = f'{bn}_linear_fit_params.txt'
    # with open(fn, 'w', newline='') as csvfile:
    #     writer = csv.writer(csvfile, dialect='JV')
    #     header = []
        
    return