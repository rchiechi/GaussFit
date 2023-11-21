'''Utility functions for the main parser class.'''

import os
import sys
import hashlib
import numpy as np
from scipy.stats import gmean
from tempfile import NamedTemporaryFile
from gaussfit.colors import RED, RS


def throwimportwarning(msg):
    print("\n\t\t%s> > > Error importing numpy/pandas/scipy! %s%s%s < < <%s" % (
          RED, RS, str(msg), RED, RS))
    print("Try pip3 install <module>")
    sys.exit()


def getdistances(coeff, X, Y):
    _a = [[], []]
    _b = [X, Y]
    for x in np.linspace(min(X), max(X), 500):
        _a[0].append(x)
        _a[1].append(linear(x, *coeff))
    A = np.array(_a)
    B = np.array(_b)
    distances = []
    # rotate three times to get J-pairs in order
    for p in np.rot90(B, k=3):
        _d = []
        for q in np.rot90(A, k=3):
            _d.append(np.sqrt((p[0] - q[0])**2 + (p[1] - q[1])**2))
        distances.append(min(_d))
    return distances


def signedgmean(Y):
    '''
    Return a geometric average with the
    same sign as the input data assuming
    all values have the same sign
    '''
    if len(Y[Y < 0]):
        Ym = -1 * gmean(abs(Y))
    else:
        Ym = gmean(abs(Y))
    return Ym


def gauss(x, *p):
    '''
    Return a gaussian function
    '''
    A, mu, sigma = p
    return A * np.exp(-(x - mu)**2 / (2. * sigma**2))


def linear(x, *p):
    '''
    Return a line
    '''
    m, b = p
    return (m * x) + b


def lorenz(x, *p):
    '''
    Return a lorenzian function
    '''
    A, mu, B = p
    return A / ((x - mu)**2 + B**2)


def printFN(logger, FN):
    '''
    Print Vtrans values to the command line for convinience
    '''
    for key in ('pos', 'neg'):
        try:
            logger.info("|Vtrans %s| Gauss-mean: %0.4f Standard Deviation: %f",
                        key, FN[key]['mean'], FN[key]['std'])
            logger.info("|Vtrans %s| Geometric-mean: %0.4f Standard Deviation: %f",
                        key, FN[key]['Gmean'], FN[key]['Gstd'])
        except KeyError:
            logger.warn("Something went wrong with Vtrans calculations.")


def gettmpfilename():
    with NamedTemporaryFile(delete=False) as __tmpfile:
        fn = __tmpfile.name
    return fn


def getfilechecksum(fn):
    file_hash = hashlib.md5()
    if os.path.exists(fn):
        with open(fn, 'rb') as fh:
            while chunk := fh.read(8192):
                file_hash.update(chunk)
    return file_hash.digest()

# def guessSegments(_df):
#     '''
#     Try to guess how many segments to expect from a dataframe
#     counting from 1.
#     '''
#     V_max = _df.V.max()
#     V_min = _df.V.min()
#     n_files = len(_df.index.levels[0])
#     n_zeros = _df.V.value_counts()[0] / n_files
#     n_V_max = _df.V.value_counts()[V_max] / n_files
#     n_V_min = _df.V.value_counts()[V_min] / n_files
#     print(n_files)
#     print(n_V_max)
#     print(n_V_min)
#     print(n_zeros)
#     if n_V_max == n_V_min == n_zeros - 1:
#         return 4
#     else:
#         return -1
#     # TODO: logic to handle cases other than 0 -> V_min/max -> 0 -> V_min/max -> 0
