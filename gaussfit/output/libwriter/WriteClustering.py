import os
import csv
import logging
import numpy as np
import pandas as pd

logger = logging.getLogger('output')


def WriteClustering(self):
    '''Write clustering results and clustered J/V curves to separate files.'''
    
    if not self.opts.cluster or not hasattr(self, 'cluster') or self.cluster['clusterer'] is None:
        return
    
    clusterer = self.cluster['clusterer']
    jv_curves = self.cluster['jv_curves']
    cluster_labels = self.cluster['clusters']
    n_clusters = self.cluster['n_clusters']
    
    if n_clusters == 0:
        logger.warning("No valid clusters found for output")
        return
    
    # Write clustering summary statistics
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile + "_clustering_summary.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Clustering Summary"])
        writer.writerow(["Total J/V curves:", len(jv_curves)])
        writer.writerow(["Number of clusters:", n_clusters])
        writer.writerow(["Clustering method:", "kmeans"])  # Default from clustering params
        writer.writerow(["Estimation method:", self.opts.cluster_estimation_method])
        writer.writerow(["Resolution:", self.opts.cluster_resolution])
        writer.writerow(["Feature dimensions:", f"{clusterer.n_bins[0]} x {clusterer.n_bins[1]} = {clusterer.n_bins[0] * clusterer.n_bins[1]}"])
        writer.writerow([""])
        writer.writerow(["Cluster", "Number of curves"])
        
        unique_labels = np.unique(cluster_labels)
        for label in unique_labels:
            if label >= 0:  # Skip noise (-1) if using DBSCAN
                count = np.sum(cluster_labels == label)
                writer.writerow([f"Cluster {label}", count])
    
    # Write cluster assignments
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile + "_cluster_assignments.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Curve_Index", "Cluster_Label"])
        for i, label in enumerate(cluster_labels):
            writer.writerow([i, label])
    
    # Get voltage range from the curves
    all_voltages = set()
    for voltage_array, _ in jv_curves:
        all_voltages.update(voltage_array)
    sorted_voltages = sorted(list(all_voltages))
    
    # Write J/V data for each cluster
    unique_labels = np.unique(cluster_labels[cluster_labels >= 0])
    
    for cluster_id in unique_labels:
        # Get curves belonging to this cluster
        cluster_mask = cluster_labels == cluster_id
        cluster_curves = [jv_curves[i] for i in range(len(jv_curves)) if cluster_mask[i]]
        
        if len(cluster_curves) == 0:
            continue
            
        # Write clustered J/V data as separate file
        _fn = os.path.join(self.opts.cluster_dir, f"{self.opts.outfile}_cluster_{cluster_id}_data.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            
            # Create header: Voltage, then J columns for each curve in cluster
            header = ["Voltage_V"]
            for i in range(len(cluster_curves)):
                header.append(f"J_curve_{i+1}_A")
            writer.writerow(header)
            
            # Write data row by row for each voltage
            for voltage in sorted_voltages:
                row = [f"{voltage:0.4f}"]
                
                for voltage_array, current_array in cluster_curves:
                    # Find current value at this voltage for this curve
                    voltage_indices = np.where(np.isclose(voltage_array, voltage, atol=1e-6))[0]
                    if len(voltage_indices) > 0:
                        current_val = current_array[voltage_indices[0]]
                        row.append(f"{current_val:0.6E}")
                    else:
                        row.append("NaN")  # No data at this voltage for this curve
                
                writer.writerow(row)
        
        # Write cluster statistics
        _fn = os.path.join(self.opts.cluster_dir, f"{self.opts.outfile}_cluster_{cluster_id}_stats.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            writer.writerow([f"Cluster {cluster_id} Statistics"])
            writer.writerow(["Number of curves:", len(cluster_curves)])
            writer.writerow([""])
            
            # Calculate statistics for each voltage point
            writer.writerow(["Voltage_V", "Mean_Current_A", "Std_Current_A", "Min_Current_A", "Max_Current_A"])
            
            for voltage in sorted_voltages:
                currents_at_voltage = []
                
                for voltage_array, current_array in cluster_curves:
                    voltage_indices = np.where(np.isclose(voltage_array, voltage, atol=1e-6))[0]
                    if len(voltage_indices) > 0:
                        currents_at_voltage.append(current_array[voltage_indices[0]])
                
                if len(currents_at_voltage) > 0:
                    currents_array = np.array(currents_at_voltage)
                    writer.writerow([
                        f"{voltage:0.4f}",
                        f"{np.mean(currents_array):0.6E}",
                        f"{np.std(currents_array):0.6E}",
                        f"{np.min(currents_array):0.6E}",
                        f"{np.max(currents_array):0.6E}"
                    ])
    
    # Save clustering plots to the clustering directory
    try:
        import matplotlib.pyplot as plt
        
        # Save main clustering results plot
        clustering_fig = clusterer.plot_results(jv_curves, figsize=(15, 10), log_scale=False)
        clustering_fig.savefig(os.path.join(self.opts.cluster_dir, f"{self.opts.outfile}_clustering_results.png"), 
                              format="png", dpi=300, bbox_inches='tight')
        plt.close(clustering_fig)
        
        # Save 2D histograms if there aren't too many clusters
        if n_clusters <= 5:
            histogram_fig = clusterer.plot_2d_histograms(jv_curves, n_examples=3)
            histogram_fig.savefig(os.path.join(self.opts.cluster_dir, f"{self.opts.outfile}_cluster_histograms.png"), 
                                 format="png", dpi=300, bbox_inches='tight')
            plt.close(histogram_fig)
    
    except Exception as e:
        logger.warning(f"Failed to save clustering plots: {e}")
    
    # Write EGaIn-compatible trace files for each cluster (for further GaussFit analysis)
    _write_egain_traces_by_cluster(self, clusterer, jv_curves, cluster_labels, n_clusters)
    
    logger.info(f"Wrote clustering results for {n_clusters} clusters to {self.opts.cluster_dir}")


def _write_egain_traces_by_cluster(writer, clusterer, jv_curves, cluster_labels, n_clusters):
    """
    Write EGaIn-compatible trace files for each cluster that can be parsed back into GaussFit.
    This reconstructs the original trace structure from the parser's data.
    """
    if not hasattr(writer.parser, 'df') or not hasattr(writer.parser, 'trace_mapping'):
        logger.warning("Original trace data or trace mapping not available for EGaIn trace reconstruction")
        return
    
    try:
        # Get the original data and trace mapping
        original_df = writer.parser.df
        trace_mapping = writer.parser.trace_mapping
        
        # Check if we're in EGaIn clustering mode
        egain_mode = hasattr(writer.parser, 'cluster') and writer.parser.cluster.get('egain_mode', False)
        
        # Get unique cluster labels (excluding noise)
        unique_labels = np.unique(cluster_labels[cluster_labels >= 0])
        
        for cluster_id in unique_labels:
            # Get trace indices that belong to this cluster
            cluster_mask = cluster_labels == cluster_id
            cluster_trace_indices = np.where(cluster_mask)[0]
            
            if len(cluster_trace_indices) == 0:
                continue
            
            # Collect all EGaIn trace data for this cluster
            cluster_traces_data = []
            
            for trace_idx in cluster_trace_indices:
                if egain_mode:
                    # In EGaIn mode, trace_idx directly maps to compressed trace index
                    compressed_trace_idx = trace_idx
                else:
                    # In sweep mode, convert individual sweep index to compressed trace index
                    # XY structure has individual sweeps (forward/reverse), trace_mapping has compressed traces
                    # If we have 10 individual sweeps from 5 EGaIn traces, mapping is: 
                    # sweep 0,1 -> trace 0; sweep 2,3 -> trace 1; etc.
                    compressed_trace_idx = trace_idx // 2  # Integer division to map pairs to single trace
                
                if compressed_trace_idx < len(trace_mapping):
                    trace_start, trace_end = trace_mapping[compressed_trace_idx]
                    
                    # Extract the original EGaIn trace data
                    try:
                        trace_data = original_df[trace_start:trace_end].copy()
                        # Only add if we haven't already added this trace (avoid duplicates)
                        if not any(trace_data.equals(existing_data) for existing_data in cluster_traces_data):
                            cluster_traces_data.append(trace_data)
                    except (KeyError, IndexError) as e:
                        logger.warning(f"Could not extract trace {trace_idx} (compressed {compressed_trace_idx}) for cluster {cluster_id}: {e}")
                        continue
                else:
                    logger.warning(f"Compressed trace index {compressed_trace_idx} (from trace {trace_idx}) >= trace_mapping length {len(trace_mapping)}, skipping")
            
            if cluster_traces_data:
                # Combine all traces for this cluster
                cluster_df = pd.concat(cluster_traces_data)
                
                # Write to EGaIn-compatible format
                output_file = os.path.join(writer.opts.cluster_dir, 
                                         f"{writer.opts.outfile}_cluster_{cluster_id}_egain_traces.txt")
                
                # Write with same format as original input files
                with open(output_file, 'w', newline='') as csvfile:
                    writer_csv = csv.writer(csvfile, dialect='JV')
                    
                    # Write header
                    writer_csv.writerow(["V", "J"])  # Standard EGaIn format
                    
                    # Write all trace data
                    for _, row in cluster_df.iterrows():
                        writer_csv.writerow([f"{row['V']:0.4f}", f"{row['J']:0.6E}"])
                
                logger.info(f"Wrote {len(cluster_trace_indices)} EGaIn traces for cluster {cluster_id} to egain_traces.txt")
    
    except Exception as e:
        logger.warning(f"Failed to write EGaIn-compatible traces: {e}")
        logger.debug("This is normal if clustering was run on non-EGaIn data")