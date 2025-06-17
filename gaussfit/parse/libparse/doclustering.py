from clustering.jvclustering import cluster_jv_curves
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