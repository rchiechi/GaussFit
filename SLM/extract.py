import sys
import logging
import csv
import numpy as np
import scipy.interpolate
from scipy.stats import linregress

csv.register_dialect('JV', delimiter='\t', quoting=csv.QUOTE_MINIMAL)


def findG(V, J, **kwargs):
    '''
    Take input J/V data and find conductance
    '''
    logger = kwargs.get('logger', getLogger())
    G = {'slope': 0, 'R': 0}
    if len(V) != len(J):
        logger.error("J/V data differ in length")
        return G
    if len(V) < 5:
        logger.error("J/V data are too short to find conductance.")
    if kwargs.get('unlog', False):
        J = np.power(10, np.array(J))
        for _i in range(0, len(V)):
            if V[_i] < 0:
                J[_i] = -1 * J[_i]
    if findMiddle(V)[1] is None:
        _n, _m = findMiddle(V)[0] - 2, findMiddle(V)[0] + 3
    else:
        _n, _m = findMiddle(V)[0] - 1, findMiddle(V)[1] + 2
    _x, _y = [], []
    for _i in range(_n, _m):
        if V[_i] == 0 and J[_i] != 0:
            continue
        _x.append(V[_i])
        _y.append(J[_i])
    try:
        _fit = linregress(_x, _y)
    except ValueError:
        logger.warning("Cannot extract conductance.")
        return G
    logger.debug(f"G:{_fit.slope:.2E} (R={_fit.rvalue ** 2:.2f})")
    G['slope'] = _fit.slope
    G['R'] = _fit.rvalue

    if kwargs.get('plot', False):
        lin = lambda x: _fit.slope * x + _fit.intercept
        import matplotlib.pyplot as plt
        xlin = np.linspace(min(_x), max(_x), 20)
        plt.plot(V, J, 'ro', ms=5)
        plt.plot(_x, _y, 'bx', ms=5)
        plt.plot(xlin, lin(xlin), 'g', lw=3)
        plt.show()

    return G


def findMiddle(input_list):
    middle = float(len(input_list)) / 2
    if middle % 2 != 0:
        return int(middle - .5), None
    else:
        return int(middle), int(middle - 1)


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
            _j = np.power(10, _j)
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

    logger.info("* * * * * * Finding Vtrans * * * * * * * *")
    FN['err'] = False
    try:
        splpos = scipy.interpolate.InterpolatedUnivariateSpline(FN['pos'][0], FN['pos'][1], k=4)
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
        if not FN['err']:
            plt.plot(neg_min_x, splneg(neg_min_x), 'bx', ms=10)
            plt.plot(pos_min_x, splpos(pos_min_x), 'bx', ms=10)
        plt.show()

    return FN


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
    G = findG(x, y, unlog=True, logger=logger, plot=True)