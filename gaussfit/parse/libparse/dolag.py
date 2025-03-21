import logging
from logging.handlers import QueueHandler
from gaussfit.parse.libparse.util import throwimportwarning, getdistances
import numpy as np
from scipy.stats import linregress




def doLag(conn, opts, que, xy):
    '''
    Make a lag plot of Y
    '''
    lag = {}
    exclude_warnings = []
    logger = logging.getLogger(__package__+".dolag")
    logger.addHandler(QueueHandler(que))
    for x, group in xy:
        lag[x] = {'lagplot': np.array([[], []]), 'filtered': np.array([])}
        Y = group['logJ']
        _lag = [[], []]
        _filtered = []
        for i in range(0, (len(Y)-len(Y) % 2), 2):
            _lag[0].append(Y.iloc[i])
            _lag[1].append(Y.iloc[i+1])
        try:
            m, b, r, _, _ = linregress(_lag[0], _lag[1])
        except FloatingPointError:
            logger.warn("Error computing lag for J = %s", _lag[0])
            continue
        # logger.debug("R-squared: %s" % (r**2))
        distances = getdistances((m, b), _lag[0], _lag[1])
        min_distance = min(distances)
        logger.debug("Distance from lag: %s", min_distance)
        if min_distance > opts.lagcutoff:
            logger.debug("Found a high degree of scatter in lag plot (%0.4f)", min_distance)
        if r**2 < 0.9:
            logger.debug("Poor line-fit to lag plot (R-squared: %s)", (r**2))
        tossed = 0
        for _t in enumerate(distances):
            if _t[1] < opts.lagcutoff:
                _filtered.append(_lag[0][_t[0]])
                _filtered.append(_lag[1][_t[0]])
            else:
                tossed += 1
        if not _filtered:
            # logger.warning("Lag filter excluded all data at %s V", x)
            exclude_warnings.append(str(x))
        if tossed > 0:
            logger.info("Lag filtered excluded %s data points at %s V", tossed, x)
        lag[x]['lagplot'] = np.array(_lag)
        lag[x]['filtered'] = np.array(_filtered)
    if exclude_warnings:
        logger.warning("Lag filter excluded all data at these voltages: %s", ",".join(exclude_warnings))
    logger.info("Lag done.")
    conn.put(lag)
