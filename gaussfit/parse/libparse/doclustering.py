from clustering.jvclustering import cluster_jv_curves
from clustering.egain_clustering import cluster_egain_traces
from gaussfit.parse.libparse.util import throwimportwarning
import logging
from logging.handlers import QueueHandler

try:
    import numpy as np
except ImportError as msg:
    throwimportwarning(msg)

def doclustering(conn, opts, que, complete_egain_traces, XY):
    try:
        conn.put(_doclustering(opts, que, complete_egain_traces, XY))
    except Exception as e:
        conn.put(e)

def _doclustering(opts, que, complete_egain_traces, XY):
    '''
    Perform clustering analysis on J/V curves
    '''
    
    if not opts.cluster:
        return 0
    logger = logging.getLogger(__package__+".clustering")
    logger.addHandler(QueueHandler(que))
    logger.info("* * * * * * Computing Clustering  * * * * * * * * *")
    # Check if we should use EGaIn-specific clustering
    if opts.cluster_as_egain:
        return _doclustering_egain(logger, opts, complete_egain_traces)
    else:
        return _doclustering_sweeps(logger, opts, XY)


def _doclustering_egain(logger, opts, egain_traces):
    '''
    Perform clustering analysis on complete EGaIn traces (0→min→0→max→0)
    '''
    cluster = {}
    logger.info("Using EGaIn trace clustering mode")
    
    if len(egain_traces) == 0:
        logger.warning("No EGaIn traces found for clustering analysis.")
        return cluster
    
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
        'resolution': opts.cluster_resolution,
        'feature_type': 'egain_complete',  # Use EGaIn-specific features
        'cluster_method': 'kmeans',
        'n_clusters': None,  # Auto-determine
        'estimation_method': opts.cluster_estimation_method,
        'dim_reduction': 'umap'
    }
    
    # Perform EGaIn clustering
    logger.info(f"Clustering {len(egain_traces)} complete EGaIn traces with resolution {opts.cluster_resolution}...")
    clusterer = cluster_egain_traces(egain_traces, clustering_params, logger=logger)
    
    # Store results in cluster (compatible with existing output system)
    cluster['clusterer'] = clusterer
    cluster['clusters'] = clusterer.labels
    cluster['jv_curves'] = egain_traces  # Store EGaIn traces as jv_curves for compatibility
    cluster['n_clusters'] = len(np.unique(clusterer.labels[clusterer.labels >= 0]))
    cluster['egain_mode'] = True  # Flag for output writers
    
    return cluster



def _doclustering_sweeps(logger, opts, XY):
    '''
    Perform clustering analysis on individual voltage sweeps (original behavior)
    '''
    logger.info("Using individual sweep clustering mode")
    cluster = {}
    # Prepare J/V curve data for clustering
    jv_curves = []
    voltage_range = None
    current_range = None
    
    # The XY structure has voltage as keys and each voltage has lists of current values
    # We need to reconstruct individual J/V traces from this structure
    voltages = sorted(list(XY.keys()))
    
    # Determine the number of traces (curves) from the first voltage point
    if voltages:
        n_traces = len(XY[voltages[0]]['Y'])
    else:
        logger.warning("No voltage data found for clustering analysis.")
        return cluster
    
    # Reconstruct individual J/V traces
    for trace_idx in range(n_traces):
        voltage_vals = []
        current_vals = []
        
        for voltage in voltages:
            if trace_idx < len(XY[voltage]['Y']):
                voltage_vals.append(voltage)
                current_vals.append(XY[voltage]['Y'][trace_idx])
        
        if len(voltage_vals) > 0:
            jv_curves.append((np.array(voltage_vals), np.array(current_vals)))
    
    if len(jv_curves) == 0:
        logger.warning("No J/V curves found for clustering analysis.")
        return cluster
    
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
        'resolution': opts.cluster_resolution,  # Use command line resolution option
        'feature_type': 'flatten',
        'cluster_method': 'kmeans',
        'n_clusters': None,  # Auto-determine
        'estimation_method': opts.cluster_estimation_method,  # Use command line option
        'dim_reduction': 'umap',
        'log_current': False  # Keep linear for consistency with existing data
    }
    
    # Perform clustering
    logger.info(f"Clustering {len(jv_curves)} J/V curves with resolution {opts.cluster_resolution}...")
    clusterer = cluster_jv_curves(jv_curves, clustering_params, logger=logger)
    
    # Store results in cluster
    cluster['clusterer'] = clusterer
    cluster['clusters'] = clusterer.labels
    cluster['jv_curves'] = jv_curves
    cluster['n_clusters'] = len(np.unique(clusterer.labels[clusterer.labels >= 0]))
    
    return cluster