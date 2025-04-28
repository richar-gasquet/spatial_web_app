#!/usr/bin/env python
"""
Batch job for clustering
"""
import os
import sys
import logging
import numpy as np
import scanpy as sc
import anndata as ad
from datetime import datetime

# Add common module to path
sys.path.append('/app/common')
from common.batch_utils import (
    update_job_status, 
    get_input_data, 
    save_output_data, 
    get_job_params
)

# Configure logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('clustering')

def process_clustering(adata, params):
    """
    Process clustering
    
    Parameters:
    -----------
    adata : anndata.AnnData
        Input data
    params : dict
        Job parameters
        
    Returns:
    --------
    anndata.AnnData
        Processed data
    dict
        Results metrics
    """
    logger.info(f"Starting clustering with parameters: {params}")
    
    # Get parameters
    resolution = params.get('resolution', 0.5)
    cluster_method = params.get('cluster_method', 'leiden')
    key_added = params.get('key_added', 'clusters')
    
    # Check if we have neighbors
    if 'neighbors' not in adata.uns:
        logger.error("No neighborhood graph found in data. Please run dimensionality reduction first.")
        raise ValueError("No neighborhood graph found in data")
    
    # Start timing
    start_time = datetime.now()
    
    # Run clustering
    if cluster_method.lower() == 'leiden':
        logger.info(f"Running Leiden clustering with resolution {resolution}")
        sc.tl.leiden(adata, resolution=resolution, key_added=key_added)
    elif cluster_method.lower() == 'louvain':
        logger.info(f"Running Louvain clustering with resolution {resolution}")
        sc.tl.louvain(adata, resolution=resolution, key_added=key_added)
    else:
        logger.error(f"Unknown clustering method: {cluster_method}")
        raise ValueError(f"Unknown clustering method: {cluster_method}")
    
    # Extract cluster information
    cluster_counts = adata.obs[key_added].value_counts().to_dict()
    n_clusters = len(cluster_counts)
    
    logger.info(f"Found {n_clusters} clusters")
    
    # Calculate processing time
    end_time = datetime.now()
    processing_time = (end_time - start_time).total_seconds()
    
    # Create results
    results = {
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'n_clusters': n_clusters,
        'cluster_counts': cluster_counts,
        'processing_time_seconds': processing_time
    }
    
    return adata, results

def main():
    """Main entry point for the clustering batch job"""
    try:
        # Update job status to running
        update_job_status('RUNNING')
        
        # Get input data
        adata = get_input_data()
        
        # Get job parameters
        params = get_job_params()
        
        # Process data
        processed_adata, results = process_clustering(adata, params)
        
        # Save output data
        output_uris = save_output_data(processed_adata, {'results': results})
        
        # Update job status to succeeded
        update_job_status('SUCCEEDED', {'output_s3_uri': output_uris['output_s3_uri']})
        
        logger.info("Job completed successfully")
        return 0
    
    except Exception as e:
        logger.error(f"Job failed: {e}", exc_info=True)
        update_job_status('FAILED', {'error': str(e)})
        return 1

if __name__ == "__main__":
    sys.exit(main())