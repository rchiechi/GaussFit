import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
from scipy.ndimage import gaussian_filter
from typing import List, Tuple, Optional, Dict, Any
import warnings
import logging
warnings.filterwarnings('ignore')

class JVCurveClustering:
    """
    Unsupervised clustering for J/V curves using 2D histogram representation
    and various clustering algorithms.
    """
    
    def __init__(self, 
                 voltage_range: Tuple[float, float] = (-2.0, 2.0),
                 current_range: Tuple[float, float] = (1e-16, 1e-3),  # Changed to 1e-16
                 n_bins: Tuple[int, int] = None,
                 resolution: int = 100,
                 log_current: bool = False,  # Changed to False for linear data
                 logger=None):
        """
        Initialize the clustering system.
        
        Parameters:
        -----------
        voltage_range : tuple
            (min, max) voltage values for the region of interest
        current_range : tuple
            (min, max) current values for the region of interest
        n_bins : tuple
            (n_voltage_bins, n_current_bins) for 2D histogram. If None, will be 
            calculated from resolution parameter.
        resolution : int
            Target resolution for feature space. Creates approximately resolution^2 
            total features for 2D histograms. Default 100 gives ~10k features.
            Historical default was ~28 (784 features), modern systems can handle 200+ (40k+ features).
        log_current : bool
            Whether to use log scale for current axis (default False for linear data)
        logger : logging.Logger, optional
            Logger instance for emitting log messages. If None, will create a default logger.
        """
        self.voltage_range = voltage_range
        self.current_range = current_range
        self.resolution = resolution
        
        # Initialize logger
        if logger is None:
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger
        
        # Calculate n_bins from resolution if not explicitly provided
        if n_bins is None:
            # Use square bins by default, but can be adjusted based on voltage/current range ratio
            voltage_span = voltage_range[1] - voltage_range[0]
            if log_current:
                current_span = np.log10(current_range[1]) - np.log10(current_range[0])
            else:
                current_span = current_range[1] - current_range[0]
            
            # Adjust bin ratio based on data range ratio, but keep total bins around resolution^2
            # Limit aspect ratio to reasonable bounds to avoid extreme rectangles
            raw_aspect_ratio = voltage_span / current_span if current_span != 0 else 1.0
            aspect_ratio = np.clip(raw_aspect_ratio, 0.2, 5.0)  # Keep ratio between 1:5 and 5:1
            
            # Calculate bins to get approximately resolution^2 total bins
            voltage_bins = int(resolution * np.sqrt(aspect_ratio))
            current_bins = int(resolution / np.sqrt(aspect_ratio))
            
            # Ensure minimum reasonable bins
            voltage_bins = max(20, voltage_bins)
            current_bins = max(20, current_bins)
            
            self.n_bins = (voltage_bins, current_bins)
        else:
            self.n_bins = n_bins
            
        self.log_current = log_current
        
        # Create bin edges
        self.voltage_bins = np.linspace(voltage_range[0], voltage_range[1], self.n_bins[0])
        if log_current:
            self.current_bins = np.logspace(np.log10(current_range[0]), 
                                          np.log10(current_range[1]), 
                                          self.n_bins[1])
        else:
            self.current_bins = np.linspace(current_range[0], current_range[1], self.n_bins[1])
            
        self.features = None
        self.labels = None
        self.scaler = StandardScaler()
        
    def curve_to_2d_histogram(self, voltage: np.ndarray, current: np.ndarray,
                            smooth: bool = True, sigma: float = 1.0) -> np.ndarray:
        """
        Convert a single J/V curve to a 2D histogram representation.
        
        Parameters:
        -----------
        voltage : array
            Voltage values
        current : array
            Current values (absolute values will be taken)
        smooth : bool
            Whether to apply Gaussian smoothing
        sigma : float
            Standard deviation for Gaussian kernel
            
        Returns:
        --------
        hist2d : array
            2D histogram of shape (n_current_bins, n_voltage_bins)
        """
        # Take absolute value of current
        current = np.abs(current)
        
        # Filter data within ROI
        mask = ((voltage >= self.voltage_range[0]) & 
                (voltage <= self.voltage_range[1]) &
                (current >= self.current_range[0]) & 
                (current <= self.current_range[1]))
        
        v_filtered = voltage[mask]
        i_filtered = current[mask]
        
        # Always return consistent shape
        if len(v_filtered) == 0:
            return np.zeros((self.n_bins[1]-1, self.n_bins[0]-1))
        
        # Create 2D histogram
        hist2d, _, _ = np.histogram2d(i_filtered, v_filtered, 
                                     bins=[self.current_bins, self.voltage_bins])
        
        # Normalize
        if hist2d.sum() > 0:
            hist2d = hist2d / hist2d.sum()
        
        # Apply smoothing if requested
        if smooth:
            hist2d = gaussian_filter(hist2d, sigma=sigma)
            
        return hist2d
    
    def extract_features(self, jv_curves: List[Tuple[np.ndarray, np.ndarray]],
                        feature_type: str = 'flatten',
                        additional_features: bool = True) -> np.ndarray:
        """
        Extract features from J/V curves.
        
        Parameters:
        -----------
        jv_curves : list of tuples
            List of (voltage, current) arrays
        feature_type : str
            'flatten': flatten 2D histogram
            'statistics': extract statistical features
            'combined': both flatten and statistics
        additional_features : bool
            Whether to add curve-specific features (max current, voltage at max, etc.)
            
        Returns:
        --------
        features : array
            Feature matrix of shape (n_curves, n_features)
        """
        features_list = []
        
        # First, determine the feature length by processing one curve
        hist2d_shape = None
        
        for idx, (voltage, current) in enumerate(jv_curves):
            # Get 2D histogram
            hist2d = self.curve_to_2d_histogram(voltage, current)
            
            if hist2d_shape is None:
                hist2d_shape = hist2d.shape
            
            features = []
            
            if feature_type in ['flatten', 'combined']:
                # Flatten the 2D histogram
                features.extend(hist2d.flatten())
            
            if feature_type in ['statistics', 'combined']:
                # Extract statistical features from histogram
                features.extend([
                    hist2d.mean(),
                    hist2d.std(),
                    hist2d.max(),
                    np.percentile(hist2d, 25),
                    np.percentile(hist2d, 50),
                    np.percentile(hist2d, 75),
                    # Center of mass
                    *np.unravel_index(hist2d.argmax(), hist2d.shape)
                ])
            
            if additional_features:
                # Add curve-specific features
                current_abs = np.abs(current)
                mask = ((voltage >= self.voltage_range[0]) & 
                       (voltage <= self.voltage_range[1]) &
                       (current_abs >= self.current_range[0]) &
                       (current_abs <= self.current_range[1]))
                
                if mask.sum() > 0:
                    v_roi = voltage[mask]
                    i_roi = current_abs[mask]
                    
                    # For linear data, don't take log
                    features.extend([
                        i_roi.max(),
                        i_roi.min(),
                        v_roi[i_roi.argmax()] if len(i_roi) > 0 else 0,
                        np.mean(i_roi),
                    ])
                else:
                    features.extend([0, 0, 0, 0])
            
            features_list.append(features)
        
        # Convert to numpy array - should now have consistent lengths
        self.features = np.array(features_list)
        return self.features
    
    def reduce_dimensions(self, method: str = 'pca', n_components: int = 2) -> np.ndarray:
        """
        Reduce dimensionality of features for visualization.
        
        Parameters:
        -----------
        method : str
            'pca', 'tsne', or 'umap'
        n_components : int
            Number of components
            
        Returns:
        --------
        reduced_features : array
            Reduced feature matrix
        """
        if self.features is None:
            raise ValueError("Extract features first using extract_features()")
        
        # Scale features
        features_scaled = self.scaler.fit_transform(self.features)
        
        if method == 'pca':
            reducer = PCA(n_components=n_components)
        elif method == 'tsne':
            reducer = TSNE(n_components=n_components, random_state=42)
        elif method == 'umap':
            reducer = umap.UMAP(n_components=n_components, random_state=42)
        else:
            raise ValueError(f"Unknown method: {method}")
        
        self.reduced_features = reducer.fit_transform(features_scaled)
        return self.reduced_features
    
    def cluster(self, method: str = 'kmeans', n_clusters: int = None, 
               estimation_method: str = 'elbow', **kwargs) -> np.ndarray:
        """
        Perform clustering on the features.
        
        Parameters:
        -----------
        method : str
            'kmeans', 'dbscan', or 'hierarchical'
        n_clusters : int
            Number of clusters (for kmeans and hierarchical). If None, will be estimated.
        estimation_method : str
            Method for estimating n_clusters: 'elbow', 'silhouette', 'gap'
        **kwargs : dict
            Additional parameters for clustering algorithms
            
        Returns:
        --------
        labels : array
            Cluster labels
        """
        if self.features is None:
            raise ValueError("Extract features first using extract_features()")
        
        # Scale features
        features_scaled = self.scaler.fit_transform(self.features)
        
        if method == 'kmeans':
            if n_clusters is None:
                # Use the new estimation method
                results = self.estimate_optimal_clusters(self.features, methods=[estimation_method])
                n_clusters = results[estimation_method]['optimal_k']
            clusterer = KMeans(n_clusters=n_clusters, n_init=100, random_state=42, **kwargs)
            
        elif method == 'dbscan':
            eps = kwargs.get('eps', 0.5)
            min_samples = kwargs.get('min_samples', 5)
            clusterer = DBSCAN(eps=eps, min_samples=min_samples)
            
        elif method == 'hierarchical':
            if n_clusters is None:
                # Use the new estimation method
                results = self.estimate_optimal_clusters(self.features, methods=[estimation_method])
                n_clusters = results[estimation_method]['optimal_k']
            clusterer = AgglomerativeClustering(n_clusters=n_clusters, **kwargs)
            
        else:
            raise ValueError(f"Unknown clustering method: {method}")
        
        self.labels = clusterer.fit_predict(features_scaled)
        self.clusterer = clusterer
        
        return self.labels
    
    def plot_results(self, jv_curves: List[Tuple[np.ndarray, np.ndarray]], 
                    figsize: Tuple[int, int] = (15, 10),
                    log_scale: bool = False) -> plt.Figure:
        """
        Comprehensive visualization of clustering results.
        
        Parameters:
        -----------
        log_scale : bool
            Whether to plot current in log scale (default False for linear data)
        """
        if self.labels is None:
            raise ValueError("Perform clustering first using cluster()")
        
        fig = plt.figure(figsize=figsize)
        
        # 1. Plot 2D projection with clusters
        ax1 = plt.subplot(2, 3, 1)
        if hasattr(self, 'reduced_features'):
            scatter = ax1.scatter(self.reduced_features[:, 0], 
                                self.reduced_features[:, 1], 
                                c=self.labels, cmap='viridis', alpha=0.6)
            ax1.set_xlabel('Component 1')
            ax1.set_ylabel('Component 2')
        ax1.set_title('2D Projection of Clusters')
        plt.colorbar(scatter, ax=ax1)
        
        # 2. Plot example curves from each cluster
        unique_labels = np.unique(self.labels)
        n_examples = 3
        
        for i, label in enumerate(unique_labels[:5]):  # Show max 5 clusters
            if label == -1:  # Skip noise in DBSCAN
                continue
                
            ax = plt.subplot(2, 3, i+2)
            
            # Get curves in this cluster
            cluster_indices = np.where(self.labels == label)[0]
            example_indices = cluster_indices[:n_examples]
            
            for idx in example_indices:
                voltage, current = jv_curves[idx]
                if log_scale:
                    ax.semilogy(voltage, np.abs(current), alpha=0.7)
                else:
                    ax.plot(voltage, current, alpha=0.7)  # Keep linear for linear data
            
            ax.set_xlabel('Voltage (V)')
            ax.set_ylabel('Current (A)' if not log_scale else '|Current| (A)')
            ax.set_title(f'Cluster {label} (n={len(cluster_indices)})')
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig
    
    def plot_2d_histograms(self, jv_curves: List[Tuple[np.ndarray, np.ndarray]], 
                          n_examples: int = 3) -> plt.Figure:
        """
        Plot example 2D histograms from each cluster.
        """
        if self.labels is None:
            raise ValueError("Perform clustering first using cluster()")
        
        unique_labels = np.unique(self.labels)
        unique_labels = unique_labels[unique_labels != -1]  # Remove noise label
        n_clusters = len(unique_labels)
        
        if n_clusters == 0:
            raise ValueError("No valid clusters found")
        
        fig, axes = plt.subplots(n_clusters, n_examples, 
                                figsize=(3*n_examples, 3*n_clusters))
        
        if n_clusters == 1:
            axes = axes.reshape(1, -1)
        if n_examples == 1:
            axes = axes.reshape(-1, 1)
        
        for i, label in enumerate(unique_labels):
            cluster_indices = np.where(self.labels == label)[0]
            example_indices = cluster_indices[:min(n_examples, len(cluster_indices))]
            
            for j, idx in enumerate(example_indices):
                voltage, current = jv_curves[idx]
                hist2d = self.curve_to_2d_histogram(voltage, current)
                
                ax = axes[i, j] if n_clusters > 1 and n_examples > 1 else axes.flat[i*n_examples + j]
                
                # Adjust extent based on whether using log scale
                if self.log_current:
                    extent = [self.voltage_range[0], self.voltage_range[1],
                             np.log10(self.current_range[0]), 
                             np.log10(self.current_range[1])]
                    ylabel = 'log10(|I|)'
                else:
                    extent = [self.voltage_range[0], self.voltage_range[1],
                             self.current_range[0], self.current_range[1]]
                    ylabel = 'Current (A)'
                
                im = ax.imshow(hist2d, aspect='auto', origin='lower',
                             extent=extent, cmap='viridis')
                
                if j == 0:
                    ax.set_ylabel(f'Cluster {label}\n{ylabel}')
                if i == n_clusters - 1:
                    ax.set_xlabel('V (V)')
                
                ax.set_title(f'Curve {idx}')
        
        plt.tight_layout()
        return fig

    def estimate_optimal_clusters(self, features: np.ndarray, 
                                methods: List[str] = ['elbow', 'silhouette', 'gap'],
                                max_k: int = 10) -> Dict[str, Any]:
        """
        Estimate optimal number of clusters using multiple methods.
        
        Parameters:
        -----------
        features : array
            Feature matrix
        methods : list
            Methods to use: 'elbow', 'silhouette', 'gap', 'davies_bouldin'
        max_k : int
            Maximum number of clusters to try
            
        Returns:
        --------
        results : dict
            Optimal k for each method and visualization data
        """
        from sklearn.metrics import silhouette_score, davies_bouldin_score
        
        # Ensure we have a reasonable range for K
        max_feasible_k = min(max_k, len(features) - 1)  # Can't have more clusters than data points - 1
        K = range(2, max(3, max_feasible_k + 1))  # Ensure at least one value in range
        results = {}
        
        # Scale features once
        features_scaled = self.scaler.fit_transform(features)
        
        # If we have too few data points for clustering, just use 2 clusters
        if len(features) < 3:
            self.logger.warning(f"Only {len(features)} data points - defaulting to 2 clusters")
            results['elbow'] = {'optimal_k': 2}
            return results
        
        if 'elbow' in methods:
            # Elbow method with better detection
            inertias = []
            for k in K:
                kmeans = KMeans(n_clusters=k, n_init=10, random_state=42)
                kmeans.fit(features_scaled)
                inertias.append(kmeans.inertia_)
            
            # Multiple elbow detection strategies
            if len(inertias) > 2:
                # Method 1: Maximum second derivative
                diff1 = np.diff(inertias)
                diff2 = np.diff(diff1)
                elbow1 = np.argmax(diff2) + 2
                
                # Method 2: Knee point detection (if kneed is available)
                try:
                    from kneed import KneeLocator
                    kn = KneeLocator(list(K), inertias, curve='convex', direction='decreasing')
                    elbow2 = kn.elbow if kn.elbow else elbow1
                except ImportError:
                    # Fall back to simple elbow detection if kneed is not available
                    elbow2 = elbow1
                
                results['elbow'] = {
                    'optimal_k': elbow2,
                    'inertias': inertias,
                    'k_values': list(K)
                }
        
        if 'silhouette' in methods:
            # Silhouette coefficient (higher is better)
            silhouette_scores = []
            for k in K:
                kmeans = KMeans(n_clusters=k, n_init=10, random_state=42)
                labels = kmeans.fit_predict(features_scaled)
                score = silhouette_score(features_scaled, labels)
                silhouette_scores.append(score)
            
            optimal_k = list(K)[np.argmax(silhouette_scores)]
            results['silhouette'] = {
                'optimal_k': optimal_k,
                'scores': silhouette_scores,
                'k_values': list(K)
            }
        
        if 'davies_bouldin' in methods:
            # Davies-Bouldin index (lower is better)
            db_scores = []
            for k in K:
                kmeans = KMeans(n_clusters=k, n_init=10, random_state=42)
                labels = kmeans.fit_predict(features_scaled)
                score = davies_bouldin_score(features_scaled, labels)
                db_scores.append(score)
            
            optimal_k = list(K)[np.argmin(db_scores)]
            results['davies_bouldin'] = {
                'optimal_k': optimal_k,
                'scores': db_scores,
                'k_values': list(K)
            }
        
        if 'gap' in methods:
            # Gap statistic
            gaps, stds = self._gap_statistic(features_scaled, max_k=max_k-1)
            
            # Find optimal k using gap statistic criterion
            for i in range(len(gaps) - 1):
                if gaps[i] >= gaps[i + 1] - stds[i + 1]:
                    optimal_k = i + 1  # +1 because index starts at 0
                    break
            else:
                optimal_k = len(gaps)
            
            results['gap'] = {
                'optimal_k': optimal_k + 1,  # +1 because we start from k=1
                'gaps': gaps,
                'stds': stds,
                'k_values': list(range(1, len(gaps) + 1))
            }
        
        return results
    
    def _gap_statistic(self, data: np.ndarray, max_k: int = 10, 
                    n_refs: int = 10) -> Tuple[List[float], List[float]]:
        """
        Calculate gap statistic for determining optimal number of clusters.
        """
        gaps = []
        stds = []
        
        for k in range(1, max_k + 1):
            # Cluster original data
            kmeans = KMeans(n_clusters=k, n_init=10, random_state=42)
            kmeans.fit(data)
            
            # Calculate within-cluster dispersion
            W_k = self._calculate_W_k(data, kmeans.labels_, kmeans.cluster_centers_)
            
            # Generate reference datasets and calculate expected W_k
            W_k_refs = []
            for _ in range(n_refs):
                # Generate random data within the bounding box of original data
                random_data = np.random.uniform(
                    low=data.min(axis=0),
                    high=data.max(axis=0),
                    size=data.shape
                )
                
                kmeans_ref = KMeans(n_clusters=k, n_init=10, random_state=42)
                kmeans_ref.fit(random_data)
                W_k_ref = self._calculate_W_k(random_data, kmeans_ref.labels_, 
                                            kmeans_ref.cluster_centers_)
                W_k_refs.append(np.log(W_k_ref))
            
            # Calculate gap statistic
            gap = np.mean(W_k_refs) - np.log(W_k)
            std = np.std(W_k_refs)
            
            gaps.append(gap)
            stds.append(std)
        
        return gaps, stds
    
    def _calculate_W_k(self, data: np.ndarray, labels: np.ndarray, 
                    centers: np.ndarray) -> float:
        """
        Calculate within-cluster sum of squares.
        """
        W_k = 0
        for i in range(len(centers)):
            cluster_points = data[labels == i]
            if len(cluster_points) > 0:
                W_k += np.sum((cluster_points - centers[i])**2)
        return W_k
    
    def plot_cluster_evaluation(self, results: Dict[str, Any], 
                            figsize: Tuple[int, int] = (15, 10)) -> plt.Figure:
        """
        Visualize different cluster evaluation metrics.
        """
        n_methods = len(results)
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        axes = axes.flatten()
        
        for idx, (method, data) in enumerate(results.items()):
            if idx >= 4:
                break
                
            ax = axes[idx]
            k_values = data['k_values']
            
            if method == 'elbow':
                ax.plot(k_values, data['inertias'], 'bo-')
                ax.axvline(x=data['optimal_k'], color='r', linestyle='--', 
                        label=f'Optimal k={data["optimal_k"]}')
                ax.set_ylabel('Inertia')
                ax.set_title('Elbow Method')
                
            elif method == 'silhouette':
                ax.plot(k_values, data['scores'], 'go-')
                ax.axvline(x=data['optimal_k'], color='r', linestyle='--',
                        label=f'Optimal k={data["optimal_k"]}')
                ax.set_ylabel('Silhouette Score')
                ax.set_title('Silhouette Analysis')
                
            elif method == 'davies_bouldin':
                ax.plot(k_values, data['scores'], 'mo-')
                ax.axvline(x=data['optimal_k'], color='r', linestyle='--',
                        label=f'Optimal k={data["optimal_k"]}')
                ax.set_ylabel('Davies-Bouldin Index')
                ax.set_title('Davies-Bouldin Method')
                
            elif method == 'gap':
                ax.errorbar(k_values, data['gaps'], yerr=data['stds'], 
                        fmt='co-', capsize=5)
                ax.axvline(x=data['optimal_k'], color='r', linestyle='--',
                        label=f'Optimal k={data["optimal_k"]}')
                ax.set_ylabel('Gap Statistic')
                ax.set_title('Gap Statistic Method')
            
            ax.set_xlabel('Number of Clusters (k)')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig

# Example usage function
def cluster_jv_curves(jv_data: List[Tuple[np.ndarray, np.ndarray]], 
                     clustering_params: Dict[str, Any] = None,
                     logger=None) -> JVCurveClustering:
    """
    Main function to cluster J/V curves.
    
    Parameters:
    -----------
    jv_data : list of tuples
        List of (voltage, current) arrays
    clustering_params : dict
        Parameters for clustering
    logger : logging.Logger, optional
        Logger instance for emitting log messages
        
    Returns:
    --------
    clusterer : JVCurveClustering
        Fitted clustering object
    """
    if clustering_params is None:
        clustering_params = {}
    
    # Initialize clusterer
    clusterer = JVCurveClustering(
        voltage_range=clustering_params.get('voltage_range', (-2.0, 2.0)),
        current_range=clustering_params.get('current_range', (1e-16, 1e-3)),
        n_bins=clustering_params.get('n_bins', None),  # Let resolution parameter calculate this
        resolution=clustering_params.get('resolution', 100),  # Default resolution
        log_current=clustering_params.get('log_current', False),  # Default to False
        logger=logger
    )
    
    # Extract features
    clusterer.logger.info("Extracting features...")
    total_bins = clusterer.n_bins[0] * clusterer.n_bins[1]
    clusterer.logger.info(f"Using resolution {clusterer.resolution} -> {clusterer.n_bins[0]} x {clusterer.n_bins[1]} = {total_bins} histogram features")
    features = clusterer.extract_features(
        jv_data, 
        feature_type=clustering_params.get('feature_type', 'flatten'),  # Changed default
        additional_features=clustering_params.get('additional_features', True)
    )
    clusterer.logger.info(f"Total feature shape: {features.shape}")
    
    # Reduce dimensions for visualization
    clusterer.logger.info("Reducing dimensions...")
    reduced = clusterer.reduce_dimensions(
        method=clustering_params.get('dim_reduction', 'umap'),
        n_components=2
    )
    
    # Perform clustering
    clusterer.logger.info("Clustering...")
    labels = clusterer.cluster(
        method=clustering_params.get('cluster_method', 'kmeans'),
        n_clusters=clustering_params.get('n_clusters', None),
        estimation_method=clustering_params.get('estimation_method', 'elbow')
    )
    
    n_clusters = len(np.unique(labels[labels >= 0]))
    clusterer.logger.info(f"Found {n_clusters} clusters")
    
    return clusterer


# Example of how to use with synthetic data - updated for linear data
if __name__ == "__main__":
    # Generate synthetic J/V curves for testing
    np.random.seed(42)
    
    def generate_linear_jv_curve(curve_type='resistor', noise_level=0.1):
        """Generate synthetic linear J/V curves"""
        voltage = np.linspace(-2, 2, 200)
        
        if curve_type == 'resistor':
            # Pure resistor
            R = 10**np.random.uniform(3, 6)  # 1k to 1M ohm
            current = voltage / R
            
        elif curve_type == 'resistor_offset':
            # Resistor with offset
            R = 10**np.random.uniform(3, 6)
            offset = np.random.uniform(-1e-9, 1e-9)
            current = voltage / R + offset
            
        elif curve_type == 'nonlinear':
            # Slightly nonlinear
            R = 10**np.random.uniform(3, 6)
            alpha = np.random.uniform(0.01, 0.1)
            current = voltage / R * (1 + alpha * voltage**2)
        
        # Add noise
        current += np.random.normal(0, noise_level * 1e-9, len(current))
        
        return voltage, current
    
    # Generate dataset with different curve types
    jv_curves = []
    true_labels = []
    
    for curve_type in ['resistor', 'resistor_offset', 'nonlinear']:
        for _ in range(20):
            v, i = generate_linear_jv_curve(curve_type, noise_level=0.01)
            jv_curves.append((v, i))
            true_labels.append(curve_type)
    
    # Run clustering
    params = {
        'voltage_range': (-2.0, 2.0),
        'current_range': (1e-16, 1e-6),  # Adjusted for linear data
        'resolution': 60,  # Use resolution parameter instead of n_bins
        'feature_type': 'flatten',
        'cluster_method': 'kmeans',
        'n_clusters': 3,
        'dim_reduction': 'umap',
        'log_current': False  # Keep linear
    }
    
    clusterer = cluster_jv_curves(jv_curves, params)
    
    # Show results
    fig = clusterer.plot_results(jv_curves, log_scale=False)  # Linear plots
    plt.show()
    
    # Show 2D histograms
    fig2 = clusterer.plot_2d_histograms(jv_curves, n_examples=3)
    plt.show()