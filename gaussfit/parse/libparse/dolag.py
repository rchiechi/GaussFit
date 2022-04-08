import pickle
import logging
from logging.handlers import QueueHandler
from multiprocessing import Process
from gaussfit.parse.libparse.util import throwimportwarning, getdistances
try:
    import numpy as np
    from scipy.stats import linregress
except ImportError as msg:
    throwimportwarning(msg)


def dolag(self, conn, que, xy):
    if self.opts.nolag:
        lag = {}
        for x, group in xy:
            lag[x] = {'lagplot': np.array([[], []]), 'filtered': np.array([])}
        return lag
    return doLag(conn, que, self.opts, xy)


def dolagmultiprocess(self, conn, que, xy):
    if self.opts.nolag:
        lag = {}
        for x, group in xy:
            lag[x] = {'lagplot': np.array([[], []]), 'filtered': np.array([])}
        return lag
    return doLag(conn, que, self.opts, xy)


class doLagMultiprocess(Process):

    def __init__(self, conn, que, opts, xy):
        super().__init__()
        self.conn = conn
        self.que = que
        self.opts = opts
        self.xy = xy

    def run(self):
        return _dolag(self.conn, self.que, self.opts, self.xy)


class doLag(doLagMultiprocess):

    # def __init__(self, conn, que, opts, xy):
    #     self.conn = conn
    #     self.que = que
    #     self.opts = opts
    #     self.xy = xy

    def start(self):
        return _dolag(self.conn, self.que, self.opts, self.xy)

    def join(self):
        return



def _dolag(conn, que, opts, xy):
    '''
    Make a lag plot of Y
    '''

    __sendattr = getattr(conn, "send", None)
    use_pipe = callable(__sendattr)
    __sendattr = getattr(conn, "write", None)
    use_pickle = callable(__sendattr)
    lag = {}

    logger = logging.getLogger(__package__+".dolag")
    logger.addHandler(QueueHandler(que))
    logger.info("* * * * * * Computing Lag  * * * * * * * * *")
    # self.loghandler.flush()
    for x, group in xy:
        lag[x] = {'lagplot': np.array([[], []]), 'filtered': np.array([])}
        Y = group['logJ']
        _lag = [[], []]
        _filtered = []
        for i in range(0, (len(Y)-len(Y) % 2), 2):
            _lag[0].append(Y[i])
            _lag[1].append(Y[i+1])
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
            logger.warning("Lag filter excluded all data at %s V", x)
        if tossed > 0:
            logger.info("Lag filtered excluded %s data points at %s V", tossed, x)
        lag[x]['lagplot'] = np.array(_lag)
        lag[x]['filtered'] = np.array(_filtered)
    logger.info("Lag done.")
    # self.loghandler.flush()
    if use_pipe:
        conn.send(lag)
        conn.close()
    else:
        with open(conn, 'w+b') as fh:
            pickle.dump(lag, fh)
        # conn.seek(0)
    # else:
    #     return lag
