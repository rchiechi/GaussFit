from gaussfit.parse.libparse.util import throwimportwarning
import logging
from logging.handlers import QueueHandler
from scipy.stats import linregress, gmean
import numpy as np

def doconductance(conn, opts, que, ohmic, avg):
    try:
        conn.put(_doconductance(opts, que, ohmic, avg))
    except Exception as e:
        conn.put(e)

def _doconductance(opts, que, ohmic, avg):
    '''
    Find the conductance using a linear regression on the first four data points.
    '''
    logger = logging.getLogger(__package__+".doconductance")
    logger.addHandler(QueueHandler(que))
    logger.info("* * * * * * Computing G * * * * * * * *")
    SLM = {'G':{}}
    tossed = 0
    voltages = []
    # Make a list of V = 0 the three values of V above and below
    for _i in range(-3, 4):
        _idx = avg.loc[0].index.tolist().index(0) + _i
        voltages.append(avg.loc[0].index.tolist()[_idx])
    for trace in avg.index.levels[0]:
        SLM['G'][trace] = np.nan
        if opts.skipohmic and trace in ohmic:
            tossed += 1
            continue
        _Y = []
        for _v in voltages:
            try:
                _Y.append(avg.loc[trace]['J'][_v])
            except KeyError:
                continue
        try:
            _fit = linregress(voltages, _Y)
        except ValueError:
            logger.warning("Cannot compute conductance (probably because of unequal voltage steps.)")
            continue
        logger.debug(f"G:{_fit.slope:.2E} (R={_fit.rvalue:.2f})")
        if _fit.rvalue ** 2 < opts.minr:
            logger.warn("Tossing G-value with R < %s", opts.minr)
            continue
        if _fit.slope > opts.maxG:
            logger.warn(f"Tossing ridiculous G-value: {_fit.slope:.2E} > {opts.maxG}")
        SLM['G'][trace] = _fit.slope
    Gavg =gmean( [SLM['G'][trace] for trace in SLM['G']], nan_policy='omit' )
    # Gavg = gmean(list(SLM['G'].values()), nan_policy='omit')
    logger.info("Average conductance: %.2E", Gavg)
    SLM['Gavg'] = Gavg
    return SLM
