import pandas as pd
import numpy as np
from collections import OrderedDict
from gaussfit.parse.libparse.util import signedgmean
from gaussfit.args import Opts

def findtraces(self):
    '''Try to find individual J/V traces. A trace is defined
    by a complete forward and reverse trace unless the input
    dataset comprises only forward or only reverse traces.
    findTraces will average the two sweeps, which may cause
    problems with very hysteretic data.'''
    def __checktraces(traces):
        if not traces:
            self.logger.error("No traces!?")  # This should never happen!
            return False
        v = self.df.loc[self.df.index.levels[0][0]]['V'].values
        st, ed = [v[0]], [v[-1]]
        for r in self.df.index.levels[0][1:]:
            try:
                s, e = self.df.loc[r]['V'].values[0], self.df.loc[r]['V'].values[-1]
            except KeyError:
                self.logger.error("Error parsing %s. Check that J/V curves are not trucated or remove file and parse again.", r)
                return False
            if s not in st or e not in ed:
                self.logger.warning(
                    'file "%s" starts and ends with weird voltages (%s -> %s)', r, s, e)
            st.append(s)
            ed.append(e)
        if Opts.tracebyfile:
            return True
        self.logger.debug("Checking V starting from slice %s:%s",
                          traces[0][0], traces[0][1])
        lt = len(self.df.V[traces[0][0]: traces[0][1]])
        for trace in traces:
            if lt != len(self.df.V[trace[0]:trace[1]]):
                # TODO Figure out how to mask/delete these files from parsing
                if trace[0][0] != trace[1][0]:
                    self.logger.warning(
                        'Unequal voltage steps somewhere bewteen "%s" (and) "%s"',
                        trace[0][0], trace[1][0])
                else:
                    self.logger.warning(
                        'Unequal voltage steps somewhere in "%s"', trace[0][0])
                self.loghandler.flush()
                return False
            self.logger.debug("Trace: %s -> %s",
                              self.df.V[trace[0]], self.df.V[trace[-1]])
        self.logger.info("Traces look good.")
        return True
    self.loghandler.flush()
    traces = []
    ntraces = 0

    if Opts.tracebyfile:
        self.logger.info("Assuming each file contains one (foward/backward) trace.")
        for r in self.df.index.levels[0]:
            traces.append(((r, self.df.loc[r].index[0]), (r, self.df.loc[r].index[-1])))
        ntraces = len(traces)

    if not ntraces and self.df.V.value_counts().index[0] == 0.0:
        # NOTE t is a tuple with both indices 0 = filename, 1 = index
        try:
            ntraces = int(self.df.V.value_counts()[0]/3)  # Three zeros in every trace!
            self.logger.info("This looks like an EGaIn dataset.")
            for t in zip(*(iter(self.df[self.df.V == 0.00].V.index),) * 3):
                # Make an iterator that aggregates elements from each of the iterables.
                # Zip is an iterator of tuples, where the i-th tuple contains the i-th element
                # from each of the argument sequences or iterables
                traces.append((t[0], t[2]))

        except ValueError as msg:
            self.logger.warning("Did not find three-zero (EGaIn) traces!")
            self.logger.debug(str(msg))
    if not ntraces:
        self.logger.warning("This does not look like an EGaIn dataset.")
        try:
            ntraces = int(self.df.V.value_counts()[self.df.V[0]]/2)  # Two end-pointss in every trace!
            for t in zip(*(iter(self.df[self.df.V == self.df.V[0]].V.index),) * 2):
                traces.append((t[0], t[1]))
        except ValueError:
            self.logger.warning("Did not find three-zero (EGaIn) traces!")

    if not ntraces or not __checktraces(traces):
        # traces contains indices where traces (start, stop) in self.df.
        self.logger.warning("Recomputing traces based on repeat values")
        traces = []
        Vinit = self.df.V[0]
        trace = []
        for row in self.df[0:].iterrows():
            if row[1].V == Vinit:
                if len(trace) <= 1:
                    trace.append(row[0])
                elif len(trace) == 2:
                    traces.append(trace)
                    trace = [row[0]]
        ntraces = len(traces)
    if not __checktraces(traces):
        self.logger.error("Problem with traces: FN and derivative probably will not work correctly!")
        self.loghandler.flush()
    self.logger.info("Found %s traces (%s).", ntraces, len(traces))
    idx = []
    frames = {}
    self.logger.info("Compressing forward/reverse sweeps to single traces.")
    for col, _t in enumerate(traces):
        fbtrace = self.df[traces[col][0]:traces[col][1]].sort_values('V')
        avg = OrderedDict({'J': [], 'FN': []})
        idx = []
        for x, group in fbtrace.groupby('V'):
            idx.append(x)
            avg['J'].append(signedgmean(group['J']))
            fn = np.mean(group['FN'])
            # fn = self.signedgmean(group['FN'])
            avg['FN'].append(fn)
        frames[col] = pd.DataFrame(avg, index=idx)
    try:
        self.avg = pd.concat(frames)
    except ValueError:
        self.error = True
        self.avg = pd.DataFrame()
        self.logger.error('Unable to parse traces.')
    if ntraces == 1:
        self.error = True
        self.logger.warning('Only parsed one trace!')
