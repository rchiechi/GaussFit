from clustering.jvclustering import cluster_jv_curves
from clustering.egain_clustering import cluster_egain_traces, extract_egain_traces_from_avg
from gaussfit.parse.libparse.util import throwimportwarning

try:
    import numpy as np
except ImportError as msg:
    throwimportwarning(msg)


def doclustering(self):
    '''
    Perform clustering analysis on J/V curves
    '''
    
    if not self.opts.cluster:
        return 0
    
    # Check if we should use EGaIn-specific clustering
    if self.opts.cluster_as_egain:
        return _doclustering_egain(self)
    else:
        return _doclustering_sweeps(self)


def _doclustering_egain(self):
    '''
    Perform clustering analysis on complete EGaIn traces (0→min→0→max→0)
    '''
    self.logger.info("Using EGaIn trace clustering mode")
    
    # Extract complete EGaIn traces using the same logic as findtraces.py
    # but WITHOUT the averaging step that creates self.avg
    if not hasattr(self, 'df') or self.df.empty:
        self.logger.warning("No raw trace data found for EGaIn clustering analysis.")
        return 0
    
    if not hasattr(self, 'trace_mapping') or not self.trace_mapping:
        self.logger.warning("No trace mapping found for EGaIn clustering analysis.")
        return 0
    
    # Use the trace boundaries found by findtraces.py to extract complete EGaIn traces
    egain_traces = _extract_complete_egain_traces_from_trace_mapping(self.df, self.trace_mapping, logger=self.logger)
    
    if len(egain_traces) == 0:
        self.logger.warning("No EGaIn traces found for clustering analysis.")
        return 0
    
    # Determine voltage and current ranges from EGaIn traces
    all_voltages = []
    all_currents = []
    for v_vals, i_vals in egain_traces:
        all_voltages.extend(v_vals)
        all_currents.extend(i_vals)
    
    voltage_range = (min(all_voltages), max(all_voltages))
    current_range = (min(np.abs(all_currents)), max(np.abs(all_currents)))
    
    # Set clustering parameters for EGaIn traces
    clustering_params = {
        'voltage_range': voltage_range,
        'current_range': current_range,
        'resolution': self.opts.cluster_resolution,
        'feature_type': 'egain_complete',  # Use EGaIn-specific features
        'cluster_method': 'kmeans',
        'n_clusters': None,  # Auto-determine
        'estimation_method': self.opts.cluster_estimation_method,
        'dim_reduction': 'umap'
    }
    
    # Perform EGaIn clustering
    self.logger.info(f"Clustering {len(egain_traces)} complete EGaIn traces with resolution {self.opts.cluster_resolution}...")
    clusterer = cluster_egain_traces(egain_traces, clustering_params, logger=self.logger)
    
    # Store results in self.cluster (compatible with existing output system)
    self.cluster['clusterer'] = clusterer
    self.cluster['clusters'] = clusterer.labels
    self.cluster['jv_curves'] = egain_traces  # Store EGaIn traces as jv_curves for compatibility
    self.cluster['n_clusters'] = len(np.unique(clusterer.labels[clusterer.labels >= 0]))
    self.cluster['egain_mode'] = True  # Flag for output writers
    
    return self.cluster['n_clusters']


def _extract_complete_egain_traces_from_trace_mapping(df, trace_mapping, logger=None):
    """
    Extract complete EGaIn traces (0→min→0→max→0) from raw data using trace mapping.
    This replicates the trace finding logic from findtraces.py but without averaging.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Raw data with complete EGaIn traces
    trace_mapping : list
        List of (start_idx, end_idx) tuples for each trace (from findtraces.py)
    logger : logging.Logger, optional
        Logger for debug messages
        
    Returns:
    --------
    egain_traces : list of tuples
        List of (voltage_array, current_array) for complete EGaIn traces
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    egain_traces = []
    
    for trace_idx, (start_idx, end_idx) in enumerate(trace_mapping):
        try:
            # Extract complete EGaIn trace from raw data using the same slicing as findtraces.py
            trace_data = df[start_idx:end_idx].copy()
            
            if len(trace_data) > 0:
                voltages = trace_data['V'].values
                currents = trace_data['J'].values
                
                logger.debug(f"Complete EGaIn trace {trace_idx}: {len(voltages)} points, V range: {voltages.min():.2f} to {voltages.max():.2f}")
                logger.debug(f"Complete EGaIn trace {trace_idx} voltage sequence: {voltages[:5]}...{voltages[-5:]}")
                
                # Verify this looks like an EGaIn trace (should have exactly 3 zeros)
                zero_count = np.sum(np.abs(voltages) < 0.01)  # Count near-zero voltages
                logger.debug(f"Complete EGaIn trace {trace_idx}: {zero_count} zero-voltage points (should be 3 for EGaIn)")
                
                if zero_count >= 2:  # Should have at least 2-3 zeros for a valid EGaIn trace
                    egain_traces.append((np.array(voltages), np.array(currents)))
                else:
                    logger.warning(f"Trace {trace_idx} doesn't look like EGaIn trace (only {zero_count} zeros)")
                    
            else:
                logger.warning(f"Empty trace data for trace {trace_idx}")
                
        except (KeyError, IndexError) as e:
            logger.warning(f"Could not extract complete EGaIn trace {trace_idx}: {e}")
            continue
    
    logger.info(f"Extracted {len(egain_traces)} complete EGaIn traces using trace mapping")
    return egain_traces


def _doclustering_sweeps(self):
    '''
    Perform clustering analysis on individual voltage sweeps (original behavior)
    '''
    self.logger.info("Using individual sweep clustering mode")
    
    # Prepare J/V curve data for clustering
    jv_curves = []
    voltage_range = None
    current_range = None
    
    # The XY structure has voltage as keys and each voltage has lists of current values
    # We need to reconstruct individual J/V traces from this structure
    voltages = sorted(list(self.XY.keys()))
    
    # Determine the number of traces (curves) from the first voltage point
    if voltages:
        n_traces = len(self.XY[voltages[0]]['Y'])
    else:
        self.logger.warning("No voltage data found for clustering analysis.")
        return 0
    
    # Reconstruct individual J/V traces
    for trace_idx in range(n_traces):
        voltage_vals = []
        current_vals = []
        
        for voltage in voltages:
            if trace_idx < len(self.XY[voltage]['Y']):
                voltage_vals.append(voltage)
                current_vals.append(self.XY[voltage]['Y'][trace_idx])
        
        if len(voltage_vals) > 0:
            jv_curves.append((np.array(voltage_vals), np.array(current_vals)))
    
    if len(jv_curves) == 0:
        self.logger.warning("No J/V curves found for clustering analysis.")
        return 0
    
    # Determine voltage and current ranges
    voltage_range = (min(voltages), max(voltages))
    
    all_currents = []
    for v_vals, i_vals in jv_curves:
        all_currents.extend(i_vals)
    current_range = (min(np.abs(all_currents)), max(np.abs(all_currents)))
    
    # Set clustering parameters based on data characteristics
    clustering_params = {
        'voltage_range': voltage_range,
        'current_range': current_range,
        'resolution': self.opts.cluster_resolution,  # Use command line resolution option
        'feature_type': 'flatten',
        'cluster_method': 'kmeans',
        'n_clusters': None,  # Auto-determine
        'estimation_method': self.opts.cluster_estimation_method,  # Use command line option
        'dim_reduction': 'umap',
        'log_current': False  # Keep linear for consistency with existing data
    }
    
    # Perform clustering
    self.logger.info(f"Clustering {len(jv_curves)} J/V curves with resolution {self.opts.cluster_resolution}...")
    clusterer = cluster_jv_curves(jv_curves, clustering_params, logger=self.logger)
    
    # Store results in self.cluster
    self.cluster['clusterer'] = clusterer
    self.cluster['clusters'] = clusterer.labels
    self.cluster['jv_curves'] = jv_curves
    self.cluster['n_clusters'] = len(np.unique(clusterer.labels[clusterer.labels >= 0]))
    
    return self.cluster['n_clusters']