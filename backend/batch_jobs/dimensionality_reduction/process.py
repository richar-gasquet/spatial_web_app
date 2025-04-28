#!/usr/bin/env python
"""
Batch job for dimensionality reduction (PCA, UMAP)
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
logger = logging.getLogger('dimensionality_reduction')

def process_dim_reduction(adata, params):
    """
    Process dimensionality reduction
    
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
    logger.info(f"Starting dimensionality reduction with parameters: {params}")
    
    # Get parameters
    n_pcs = params.get('n_pcs', 50)
    n_neighbors = params.get('n_neighbors', 15)
    use_hvg = params.get('use_hvg', True)
    run_pca = params.get('run_pca', True)
    run_neighbors = params.get('run_neighbors', True)
    run_umap = params.get('run_umap', True)
    
    # Start timing
    start_time = datetime.now()
    
    # Subset to highly variable genes if requested
    if use_hvg and 'highly_variable' in adata.var.columns:
        logger.info("Using only highly variable genes")
        adata_hvg = adata[:, adata.var.highly_variable]
    else:
        logger.info("Using all genes")
        adata_hvg = adata
    
    # Run PCA
    if run_pca:
        logger.info(f"Running PCA with {n_pcs} components")
        sc.tl.pca(adata_hvg, n_comps=n_pcs, svd_solver='arpack')
        
        # Copy results back to main object if we're using a subset
        if adata_hvg is not adata:
            adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
            adata.varm['PCs'] = adata_hvg.varm['PCs']
            adata.uns['pca'] = adata_hvg.uns['pca']
    
    # Run neighbors
    if run_neighbors:
        logger.info(f"Computing neighborhood graph with {n_neighbors} neighbors")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # Run UMAP
    if run_umap and 'neighbors' in adata.uns:
        logger.info("Running UMAP")
        sc.tl.umap(adata)
    
    # Calculate processing time
    end_time = datetime.now()
    processing_time = (end_time - start_time).total_seconds()
    
    # Create results
    results = {
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'processing_time_seconds': processing_time
    }
    
    # Add dimensionality reduction info
    if 'X_pca' in adata.obsm:
        results['pca_shape'] = adata.obsm['X_pca'].shape
    
    if 'X_umap' in adata.obsm:
        results['umap_shape'] = adata.obsm['X_umap'].shape
    
    return adata, results

def main():
    """Main entry point for the dimensionality reduction batch job"""
    try:
        # Update job status to running
        update_job_status('RUNNING')
        
        # Get input data
        adata = get_input_data()
        
        # Get job parameters
        params = get_job_params()
        
        # Process data
        processed_adata, results = process_dim_reduction(adata, params)
        
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