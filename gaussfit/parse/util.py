'''Utility functions for the main parser class.'''

import numpy as np
from scipy.stats import gmean

def getdistances(coeff, X, Y):
    _a = [[],[]]
    _b = [X,Y]
    for x in np.linspace(min(X), max(X), 500):
        _a[0].append(x)
        _a[1].append(linear(x,*coeff))
    A = np.array(_a)
    B = np.array(_b)
    distances = []
    # rotate three times to get J-pairs in order
    for p in np.rot90(B,k=3):
        _d = []
        for q in np.rot90(A,k=3):
            _d.append( np.sqrt( (p[0]-q[0])**2 + (p[1]-q[1])**2 ) )
        distances.append(min(_d))
    return distances


def signedgmean(Y):
    '''
    Return a geometric average with the
    same sign as the input data assuming
    all values have the same sign
    '''
    if len(Y[Y<0]):
        Ym = -1*gmean(abs(Y))
    else: Ym = gmean(abs(Y))
    return Ym

def gauss(x, *p):
    '''
    Return a gaussian function
    '''
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def linear(x, *p):
    '''
    Return a line
    '''
    m,b = p
    return (m*x)+b

def lorenz(x, *p):
    '''
    Return a lorenzian function
    '''
    A, mu, B = p
    return A/( (x-mu)**2 + B**2)

def printFN(logger, FN):
    '''
    Print Vtrans values to the command line for convinience
    '''
    for key in ('pos', 'neg'):
        logger.info("|Vtrans %s| Gauss-mean: %0.4f Standard Deviation: %f",
            key, FN[key]['mean'], FN[key]['std'] )
        logger.info("|Vtrans %s| Geometric-mean: %0.4f Standard Deviation: %f",
            key, FN[key]['Gmean'], FN[key]['Gstd'] )
