import sys
import logging
import csv
import numpy as np
import scipy.interpolate


csv.register_dialect('JV', delimiter='\t', quoting=csv.QUOTE_MINIMAL)


def findvtrans(V, J, **kwargs):
    '''
    Take input J/V data and find Vtrans
    '''
    logger = kwargs.get('logger', getLogger())
    FN = {'err': True, 'pos': [[], []], 'neg': [[], []], 'vt_pos': 1e-15, 'vt_neg': -1e-15}
    if len(V) != len(J):
        logger.error("J/V data differ in length")
        return FN
    for _i in range(0, len(V)):
        if not V[_i]:
            continue
        _j = J[_i]
        if kwargs.get('unlog', False):
            _j = np.power(_j, 10)
        _fn = abs(V[_i] ** 2 / _j)
        if V[_i] < 0:
            FN['neg'][0].append(V[_i])
            FN['neg'][1].append(_fn)
        if V[_i] > 0:
            FN['pos'][0].append(V[_i])
            FN['pos'][1].append(_fn)
    for _i in (0, 1):
        FN['pos'][_i] = np.array(FN['pos'][_i])
        FN['neg'][_i] = np.array(FN['neg'][_i])

    logger.info("* * * * * * Computing Vtrans * * * * * * * *")
    FN['err'] = False
    try:
        splpos = scipy.interpolate.UnivariateSpline(FN['pos'][0], FN['pos'][1], k=4)
        pos_min_x = list(np.array(splpos.derivative().roots()))
        logger.debug("Found positive vals: %s", pos_min_x)
        for _x in pos_min_x:
            if splpos(_x) == max(splpos(pos_min_x)):
                FN['vt_pos'] = _x
    except Exception as msg:
        logger.warning('Error finding derivative of FN(+) %s', str(msg))
        FN['err'] = True
    try:
        splneg = scipy.interpolate.UnivariateSpline(FN['neg'][0], FN['neg'][1], k=4)
        neg_min_x = list(np.array(splneg.derivative().roots()))
        logger.debug("Found negative vals: %s", neg_min_x)
        for _x in neg_min_x:
            if splneg(_x) == max(splneg(neg_min_x)):
                FN['vt_neg'] = _x
    except Exception as msg:
        logger.warning('Error finding derivative of FN(–) %s', str(msg))
        FN['err'] = True

    logger.debug(f'Vtrans(+): {FN["vt_pos"]:0.2f} V')
    logger.debug(f'Vtrans(-): {FN["vt_neg"]:0.2f} V')

    if kwargs.get('plot', False):
        import matplotlib.pyplot as plt
        xneg, xpos = np.linspace(min(V), -1e-5, 250), np.linspace(1e-5, max(V), 250)
        plt.plot(FN['neg'][0], FN['neg'][1], 'ro', ms=5)
        plt.plot(FN['pos'][0], FN['pos'][1], 'ro', ms=5)
        plt.plot(xneg[:-10], splneg(xneg[:-10]), 'g', lw=3)
        plt.plot(xpos[10:], splpos(xpos[10:]), 'g', lw=3)
        plt.plot(neg_min_x, splneg(neg_min_x), 'bx', ms=10)
        plt.plot(pos_min_x, splpos(pos_min_x), 'bx', ms=10)
        plt.show()

    return FN

    # print(xneg)
    # pirnt(yneg)
    # try:
    #     splpos = scipy.interpolate.UnivariateSpline(xpos, ypos, k=4)
    #     pos_min_x = list(np.array(splpos.derivative().roots()))
    #     logger.info("Found positive vals: %s", pos_min_x)
    # except Exception as msg:
    #     logger.warning('Error finding derivative of FN(+) %s', str(msg))
    #     err = True
    # try:
    #     splneg = scipy.interpolate.UnivariateSpline(xneg, yneg, k=4)
    #     neg_min_x = list(np.array(splneg.derivative().roots()))
    #     logger.info("Found negative vals: %s", neg_min_x)
    # except Exception as msg:
    #     logger.warning('Error finding derivative of FN(–) %s', str(msg))
    #     err = True
    # if err:
    #     logger.error('Cannot conintue')
    #     return None, None

    # return neg_min_x[-1], pos_min_x[-1]


def getLogger(name=__name__):
    # create logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    # create formatter
    formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)

    return logger


if __name__ == '__main__':
    logger = getLogger('Find Vtrans')
    logger.info("Parsing %s", sys.argv[1])
    x, y = [], []
    with open(sys.argv[1], newline='') as fh:
        reader = csv.reader(fh, 'JV')
        for row in reader:
            try:
                x.append(float(row[0]))
                y.append(float(row[1]))
            except ValueError:
                logger.info('Skipping row %s', row)
    FN = findvtrans(x, y, unlog=True, logger=logger, plot=True)



#     logger.info('Finding minimum of interpolated FN plot.')
#     splpos = scipy.interpolate.interp1d(FN['pos'][0], FN['pos'][1],
#                                         kind='linear', fill_value='extrapolate')
#     splneg = scipy.interpolate.interp1d(FN['neg'][0], FN['neg'][1],
#                                         kind='linear', fill_value='extrapolate')
#     xy = {'X': [], 'Y': []}
#     xneg, xpos = np.linspace(vmin, 0, 250), np.linspace(vmax, 0, 250)
#     for x in xneg:
#         if not np.isfinite(x):
#             logger.warning('Bad voltage %s', x)
#             continue
#         xy['Y'].append(float(splneg(x)))
#         xy['X'].append(float(x))
#     for x in xpos:
#         if not np.isfinite(x):
#             logger.warning('Bad voltage %s', x)
#             continue
#         xy['Y'].append(float(splpos(x)))
#         xy['X'].append(float(x))
#     fndf = pd.DataFrame(xy)
#     pidx = fndf[fndf.X > 0]['Y'].idxmin()
#     nidx = fndf[fndf.X < 0]['Y'].idxmin()
#     if err[0]:
#         try:
#             splpos = scipy.interpolate.UnivariateSpline(
#                 fndf['X'][pidx-20:pidx+20].values, fndf['Y'][pidx-20:pidx+20].values, k=4)
#             pos_min_x += list(np.array(splpos.derivative().roots()))
#         except Exception as msg:  # pylint: disable=broad-except
#             # TODO: Figure out how to catch all the scipy errors
#             self.logger.warning('Error finding FN(+) minimum from interpolated derivative, falling back to minimum. %s', str(msg))
#             pos_min_x.append(np.mean(fndf['X'][pidx-20:pidx+20].values))
#     if err[1]:
#         try:
#             splneg = scipy.interpolate.UnivariateSpline(fndf['X'][nidx-20:nidx+20].values,
#                                                         fndf['Y'][nidx-20:nidx+20].values, k=4)
#             neg_min_x += list(np.array(splneg.derivative().roots()))
#         except Exception as msg:  # pylint: disable=broad-except
#             # TODO: Figure out how to catch all the scipy errors
#             self.logger.warning('Error finding FN(–) minimum from interpolated derivative, falling back to minimum. %s', str(msg))
#             neg_min_x.append(np.mean(fndf['X'][nidx-20:nidx+20].values))

# self.SLM['Vtpos'][trace] = pos_min_x[-1]
# self.SLM['Vtneg'][trace] = neg_min_x[-1]
# except (IndexError, ValueError) as msg:
# self.logger.warning("Error finding minimum in trace: %s", str(msg))

#     if tossed:
#         self.logger.warning("Tossed %d compliance traces during FN calculation.", tossed)
#     neg_min_x = np.array(list(filter(np.isfinite, neg_min_x)))
#     pos_min_x = np.array(list(filter(np.isfinite, pos_min_x)))
#     if not len(neg_min_x) or not len(pos_min_x):
#         self.logger.error("Did not parse any FN values!")
#         self.FN["neg"], self.FN["pos"] = {}, {}
#     else:
#         self.FN["neg"] = dohistogram(self.logqueue, neg_min_x, label="Vtrans(-)", density=True)
#         self.FN["pos"] = dohistogram(self.logqueue, pos_min_x, label="Vtrans(+)", density=True)
