import pandas as pd
import numpy as np
from collections import OrderedDict
from gaussfit.parse.libparse.util import signedgmean

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
        if self.opts.tracebyfile:
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

    if self.opts.tracebyfile:
        self.logger.info("Assuming each file contains one (foward/backward) trace.")
        for r in self.df.index.levels[0]:
            traces.append(((r, self.df.loc[r].index[0]), (r, self.df.loc[r].index[-1])))
        ntraces = len(traces)

    if not ntraces and self.df.V.value_counts().index[0] == 0.0:
        # NOTE t is a tuple with both indices 0 = filename, 1 = index
        try:
            if self.df.V.value_counts()[0] % 3:
                raise ValueError("Number of zeros in Voltages not divisble by 3.")
            ntraces = int(self.df.V.value_counts()[0]/3)  # Three zeros in every trace!
            self.logger.info("This looks like an EGaIn dataset.")
            
            # Robust trace finding that preserves original dataframe order
            zero_indices = self.df[self.df.V == 0.00].index.tolist()
            # Keep zeros in original dataframe order (already sorted by default)
            
            # Process zeros sequentially in groups of 3, validating each trace
            i = 0
            while i + 2 < len(zero_indices):
                start_idx = zero_indices[i]
                middle_idx = zero_indices[i + 1] 
                end_idx = zero_indices[i + 2]
                
                # Get voltage values for this potential trace
                trace_subset = self.df.loc[start_idx:end_idx]
                v_values = trace_subset['V'].values
                
                # Validate this is a complete 0->vmax->0->vmin->0 trace
                has_positive = np.any(v_values > 0.001)  # Small tolerance for floating point
                has_negative = np.any(v_values < -0.001)
                
                # Check that the three zeros are reasonably spaced (not all consecutive)
                trace_length = len(trace_subset)
                min_trace_length = 5  # Minimum points for a meaningful trace
                
                if has_positive and has_negative and trace_length >= min_trace_length:
                    # This is a valid complete trace
                    traces.append((start_idx, end_idx))
                    i += 3  # Move to next set of 3 zeros
                else:
                    # This set of 3 zeros doesn't form a valid trace
                    # Log the issue and try shifting by 1 to realign
                    filename = start_idx[0] if hasattr(start_idx, '__len__') else 'unknown'
                    if not has_positive or not has_negative:
                        self.logger.warning(f"Incomplete trace starting at {start_idx}: missing positive or negative voltages")
                    if trace_length < min_trace_length:
                        self.logger.warning(f"Trace starting at {start_idx} too short ({trace_length} points)")
                    
                    i += 1  # Try next alignment
            
            ntraces = len(traces)
            if ntraces == 0:
                raise ValueError("No complete EGaIn traces found with proper voltage progression")

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
        
        # Enhanced fallback for inconsistent voltage ranges that preserves order
        if self.df.V.value_counts().index[0] == 0.0:
            # Try EGaIn-specific fallback for inconsistent data
            self.logger.info("Attempting robust EGaIn trace detection for inconsistent data")
            zero_indices = self.df[abs(self.df.V) < 0.001].index.tolist()  # More tolerant zero detection
            # zero_indices is already in dataframe order, preserve it
            
            # Process zeros sequentially with more tolerant validation
            i = 0
            while i + 2 < len(zero_indices):
                start_idx = zero_indices[i]
                end_idx = zero_indices[i + 2]
                
                # Verify trace completeness by checking voltage range
                trace_v = self.df.loc[start_idx:end_idx]['V']
                if len(trace_v) > 5:  # Minimum trace length
                    vmax = trace_v.max()
                    vmin = trace_v.min()
                    if vmax > 0.01 and vmin < -0.01:  # Has both positive and negative
                        traces.append((start_idx, end_idx))
                        i += 3  # Move to next group of 3
                    else:
                        self.logger.debug(f"Skipping incomplete trace at {start_idx}: vmax={vmax:.3f}, vmin={vmin:.3f}")
                        i += 1  # Try next alignment
                else:
                    self.logger.debug(f"Skipping short trace at {start_idx}: length={len(trace_v)}")
                    i += 1  # Try next alignment
        
        # Original fallback method if EGaIn-specific fallback didn't work
        if not traces:
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
    # Store trace mapping for clustering reconstruction
    self.trace_mapping = traces.copy()
    
    # DEBUG
    # from rich import print
    
    # Store complete EGaIn traces if EGaIn clustering is enabled
    self.complete_egain_traces = []
    
    self.logger.info("Compressing forward/reverse sweeps to single traces.")
    for col, _t in enumerate(traces):
        trace_data = self.df[traces[col][0]:traces[col][1]].copy()
        if len(trace_data) > 0:
                voltages = trace_data['V'].values
                currents = trace_data['J'].values
                self.complete_egain_traces.append((np.array(voltages), np.array(currents)))
        fbtrace = trace_data.sort_values('V')
        # DEBUG
        # print(f"fbtrace: {fbtrace['V']}")
        avg = OrderedDict({'J': [], 'FN': []})
        idx = []
        for x, group in fbtrace.groupby('V'):
            idx.append(x)
            avg['J'].append(signedgmean(group['J']))
            fn = np.mean(group['FN'])
            # fn = self.signedgmean(group['FN'])
            avg['FN'].append(fn)
        frames[col] = pd.DataFrame(avg, index=idx)
    self.logger.info(f"Stored {len(self.complete_egain_traces)} complete EGaIn traces for clustering.")

    try:
        self.avg = pd.concat(frames)
    except ValueError:
        self.error = True
        self.avg = pd.DataFrame()
        self.logger.error('Unable to parse traces.')
    if ntraces == 1:
        self.error = True
        self.logger.warning('Only parsed one trace!')
