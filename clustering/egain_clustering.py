"""
EGaIn-specific clustering functionality.

This module provides clustering based on complete EGaIn traces (0→min→0→max→0)
rather than individual voltage sweeps. This approach captures the full electrical
behavior of molecular junctions in a single feature representation.
"""

import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Any
import logging
from .jvclustering import JVCurveClustering, cluster_jv_curves


def extract_egain_traces_from_avg(avg_data, logger=None):
    """
    Extract complete EGaIn traces from the compressed average data structure.
    
    Parameters:
    -----------
    avg_data : pandas.DataFrame
        Multi-index DataFrame with trace indices and voltage values
    logger : logging.Logger, optional
        Logger for status messages
        
    Returns:
    --------
    egain_traces : list of tuples
        List of (voltage_array, current_array) for each complete EGaIn trace
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    egain_traces = []
    
    # Extract individual traces from the multi-level index structure
    for trace_idx in avg_data.index.get_level_values(0).unique():
        try:
            trace_data = avg_data.loc[trace_idx]
            voltages = trace_data.index.values  # Voltage values
            currents = trace_data['J'].values   # Current values
            
            if len(voltages) > 0 and len(currents) > 0:
                egain_traces.append((np.array(voltages), np.array(currents)))
                logger.debug(f"Extracted EGaIn trace {trace_idx}: {len(voltages)} points")
        except (KeyError, IndexError) as e:
            logger.warning(f"Could not extract EGaIn trace {trace_idx}: {e}")
            continue
    
    logger.info(f"Extracted {len(egain_traces)} complete EGaIn traces for clustering")
    return egain_traces


def extract_egain_features(voltage, current, voltage_range=(-2.0, 2.0), 
                          current_range=(1e-16, 1e-3), n_bins=(50, 50)):
    """
    Extract features specifically designed for complete EGaIn traces.
    
    Parameters:
    -----------
    voltage : array
        Voltage values for complete EGaIn trace
    current : array  
        Current values for complete EGaIn trace
    voltage_range : tuple
        Min/max voltage range for features
    current_range : tuple
        Min/max current range for features
    n_bins : tuple
        Number of bins for 2D histogram (voltage_bins, current_bins)
        
    Returns:
    --------
    features : array
        Feature vector representing the EGaIn trace characteristics
    """
    features = []
    
    # Basic 2D histogram representation
    current_abs = np.abs(current)
    
    # Filter data within ROI
    mask = ((voltage >= voltage_range[0]) & 
            (voltage <= voltage_range[1]) &
            (current_abs >= current_range[0]) & 
            (current_abs <= current_range[1]))
    
    v_filtered = voltage[mask]
    i_filtered = current_abs[mask]
    
    if len(v_filtered) > 0:
        # Create 2D histogram
        voltage_bins = np.linspace(voltage_range[0], voltage_range[1], n_bins[0])
        current_bins = np.linspace(current_range[0], current_range[1], n_bins[1])
        
        hist2d, _, _ = np.histogram2d(i_filtered, v_filtered, 
                                     bins=[current_bins, voltage_bins])
        
        # Normalize histogram
        if hist2d.sum() > 0:
            hist2d = hist2d / hist2d.sum()
            
        # Add flattened histogram as primary features
        features.extend(hist2d.flatten())
    else:
        # No valid data - use zeros
        features.extend(np.zeros(n_bins[0] * n_bins[1]))
    
    # EGaIn-specific features
    try:
        # Full trace characteristics
        voltage_span = voltage.max() - voltage.min()
        current_span = current_abs.max() - current_abs.min()
        
        # Zero-crossing behavior (important for EGaIn)
        zero_indices = np.where(np.abs(voltage) < 0.05)[0]  # Near-zero voltages
        if len(zero_indices) >= 3:
            # Should have 3 zeros in EGaIn trace
            zero_currents = current[zero_indices]
            zero_consistency = np.std(zero_currents)
        else:
            zero_consistency = 0
        
        # Forward/reverse asymmetry
        positive_v_mask = voltage > 0
        negative_v_mask = voltage < 0
        
        if np.any(positive_v_mask) and np.any(negative_v_mask):
            pos_max_current = np.max(np.abs(current[positive_v_mask]))
            neg_max_current = np.max(np.abs(current[negative_v_mask]))
            asymmetry_ratio = pos_max_current / (neg_max_current + 1e-15)
        else:
            asymmetry_ratio = 1.0
        
        # Maximum current positions
        max_current_idx = np.argmax(current_abs)
        voltage_at_max_current = voltage[max_current_idx] if len(voltage) > 0 else 0
        
        # Current range characteristics
        log_current_range = np.log10(current_span + 1e-15)
        
        # Add EGaIn-specific features
        features.extend([
            voltage_span,
            log_current_range,
            zero_consistency,
            asymmetry_ratio,
            voltage_at_max_current,
            current_abs.max(),
            current_abs.min(),
            np.mean(current_abs),
            np.std(current_abs)
        ])
        
    except Exception as e:
        # Fallback: add zero features if calculation fails
        features.extend([0] * 9)
    
    return np.array(features)


class EGaInClustering(JVCurveClustering):
    """
    Specialized clustering for complete EGaIn traces.
    
    Inherits from JVCurveClustering but uses EGaIn-specific feature extraction
    and clustering approaches designed for complete 0→min→0→max→0 traces.
    """
    
    def __init__(self, 
                 voltage_range: Tuple[float, float] = (-2.0, 2.0),
                 current_range: Tuple[float, float] = (1e-16, 1e-3),
                 n_bins: Tuple[int, int] = None,
                 resolution: int = 100,
                 logger=None):
        """
        Initialize EGaIn clustering with specialized parameters.
        
        Parameters:
        -----------
        voltage_range : tuple
            (min, max) voltage values for feature extraction
        current_range : tuple
            (min, max) current values for feature extraction  
        n_bins : tuple
            (n_voltage_bins, n_current_bins) for 2D histogram
        resolution : int
            Target resolution - will calculate n_bins if not provided
        logger : logging.Logger, optional
            Logger instance
        """
        super().__init__(
            voltage_range=voltage_range,
            current_range=current_range,
            n_bins=n_bins,
            resolution=resolution,
            log_current=False,  # Use linear current for EGaIn
            logger=logger
        )
        
        # EGaIn-specific settings
        self.egain_mode = True
        
    def extract_features(self, egain_traces: List[Tuple[np.ndarray, np.ndarray]],
                        feature_type: str = 'egain_complete') -> np.ndarray:
        """
        Extract features from complete EGaIn traces.
        
        Parameters:
        -----------
        egain_traces : list of tuples
            List of (voltage, current) arrays for complete EGaIn traces
        feature_type : str
            Type of features to extract ('egain_complete', 'histogram_only')
            
        Returns:
        --------
        features : array
            Feature matrix of shape (n_traces, n_features)
        """
        features_list = []
        
        for idx, (voltage, current) in enumerate(egain_traces):
            if feature_type == 'egain_complete':
                # Use specialized EGaIn feature extraction
                features = extract_egain_features(
                    voltage, current,
                    voltage_range=self.voltage_range,
                    current_range=self.current_range,
                    n_bins=self.n_bins
                )
            elif feature_type == 'histogram_only':
                # Use only 2D histogram (from parent class)
                hist2d = self.curve_to_2d_histogram(voltage, current)
                features = hist2d.flatten()
            else:
                raise ValueError(f"Unknown feature type: {feature_type}")
                
            features_list.append(features)
        
        # Convert to numpy array
        self.features = np.array(features_list)
        self.logger.info(f"Extracted EGaIn features: {self.features.shape}")
        
        return self.features


def cluster_egain_traces(egain_traces: List[Tuple[np.ndarray, np.ndarray]], 
                        clustering_params: Dict[str, Any] = None,
                        logger=None) -> EGaInClustering:
    """
    Main function to cluster complete EGaIn traces.
    
    Parameters:
    -----------
    egain_traces : list of tuples
        List of (voltage, current) arrays for complete EGaIn traces
    clustering_params : dict
        Parameters for clustering (similar to cluster_jv_curves)
    logger : logging.Logger, optional
        Logger instance
        
    Returns:
    --------
    clusterer : EGaInClustering
        Fitted clustering object
    """
    if clustering_params is None:
        clustering_params = {}
    
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Initialize EGaIn clusterer
    clusterer = EGaInClustering(
        voltage_range=clustering_params.get('voltage_range', (-2.0, 2.0)),
        current_range=clustering_params.get('current_range', (1e-16, 1e-3)),
        n_bins=clustering_params.get('n_bins', None),
        resolution=clustering_params.get('resolution', 100),
        logger=logger
    )
    
    # Extract EGaIn-specific features
    logger.info("Extracting EGaIn-specific features...")
    total_bins = clusterer.n_bins[0] * clusterer.n_bins[1]
    logger.info(f"Using resolution {clusterer.resolution} -> {clusterer.n_bins[0]} x {clusterer.n_bins[1]} = {total_bins} histogram features + 9 EGaIn features")
    
    features = clusterer.extract_features(
        egain_traces, 
        feature_type=clustering_params.get('feature_type', 'egain_complete')
    )
    logger.info(f"Total EGaIn feature shape: {features.shape}")
    
    # Reduce dimensions for visualization
    logger.info("Reducing dimensions...")
    reduced = clusterer.reduce_dimensions(
        method=clustering_params.get('dim_reduction', 'umap'),
        n_components=2
    )
    
    # Perform clustering
    logger.info("Clustering EGaIn traces...")
    labels = clusterer.cluster(
        method=clustering_params.get('cluster_method', 'kmeans'),
        n_clusters=clustering_params.get('n_clusters', None),
        estimation_method=clustering_params.get('estimation_method', 'elbow')
    )
    
    n_clusters = len(np.unique(labels[labels >= 0]))
    logger.info(f"Found {n_clusters} EGaIn trace clusters")
    
    return clusterer