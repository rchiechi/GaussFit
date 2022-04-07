import pickle
from gaussfit.parse.libparse.util import throwimportwarning, getdistances
try:
    import numpy as np
    from scipy.stats import linregress
except ImportError as msg:
    throwimportwarning(msg)


def dolag(self, conn, xy):
    '''
    Make a lag plot of Y
    '''

    __sendattr = getattr(conn, "send", None)
    use_pipe = callable(__sendattr)
    __sendattr = getattr(conn, "write", None)
    use_pickle = callable(__sendattr)
    lag = {}
    if self.opts.nolag:
        for x, group in xy:
            lag[x] = {'lagplot': np.array([[], []]), 'filtered': np.array([])}
        return lag
    self.logger.info("* * * * * * Computing Lag  * * * * * * * * *")
    self.loghandler.flush()
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
            self.logger.warn("Error computing lag for J = %s", _lag[0])
            continue
        # self.logger.debug("R-squared: %s" % (r**2))
        distances = getdistances((m, b), _lag[0], _lag[1])
        min_distance = min(distances)
        self.logger.debug("Distance from lag: %s", min_distance)
        if min_distance > self.opts.lagcutoff:
            self.logger.debug("Found a high degree of scatter in lag plot (%0.4f)", min_distance)
        if r**2 < 0.9:
            self.logger.debug("Poor line-fit to lag plot (R-squared: %s)", (r**2))
        tossed = 0
        for _t in enumerate(distances):
            if _t[1] < self.opts.lagcutoff:
                _filtered.append(_lag[0][_t[0]])
                _filtered.append(_lag[1][_t[0]])
            else:
                tossed += 1
        if not _filtered:
            self.logger.warning("Lag filter excluded all data at %s V", x)
        if tossed > 0:
            self.logger.info("Lag filtered excluded %s data points at %s V", tossed, x)
        lag[x]['lagplot'] = np.array(_lag)
        lag[x]['filtered'] = np.array(_filtered)
    self.logger.info("Lag done.")
    self.loghandler.flush()
    if use_pipe:
        conn.send(lag)
        conn.close()
    elif use_pickle:
        with open(conn.name, 'w+b') as fh:
            pickle.dump(lag, fh)
    else:
        return lag
