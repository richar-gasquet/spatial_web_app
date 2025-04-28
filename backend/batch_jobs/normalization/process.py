#!/usr/bin/env python
"""
Batch job for normalization and HVG selection
"""
import os
import sys
import logging
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
logger = logging.getLogger('normalization')

def process_normalization(adata, params):
    """
    Process normalization and HVG selection
    
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
    logger.info(f"Starting normalization with parameters: {params}")
    
    # Get parameters
    normalize = params.get('normalize', True)
    log_transform = params.get('log_transform', True)
    n_top_genes = params.get('n_top_genes', 2000)
    use_batch = params.get('use_batch', False)
    batch_key = params.get('batch_key', None)
    subset_hvg = params.get('subset_hvg', True)
    scale_data = params.get('scale_data', False)  # Default False for memory reasons
    
    # Start timing
    start_time = datetime.now()
    
    # Store original counts if not present
    if 'counts' not in adata.layers:
        logger.info("Storing original counts in layers['counts']")
        adata.layers['counts'] = adata.X.copy()
    
    # Normalize
    if normalize:
        logger.info("Normalizing total counts per cell")
        sc.pp.normalize_total(adata, inplace=True)
    
    # Log transform
    if log_transform:
        logger.info("Log-transforming the data")
        sc.pp.log1p(adata)
    
    # Find highly variable genes
    logger.info(f"Identifying {n_top_genes} highly variable genes")
    
    if use_batch and batch_key in adata.obs.columns:
        logger.info(f"Using batch effect correction with {batch_key}")
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat",
            n_top_genes=n_top_genes,
            batch_key=batch_key
        )
    else:
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat",
            n_top_genes=n_top_genes
        )
    
    # Save raw for later use
    logger.info("Storing raw counts")
    adata.raw = adata
    
    # Subset to highly variable genes if requested
    if subset_hvg:
        logger.info("Subsetting to highly variable genes")
        if use_batch and 'highly_variable_nbatches' in adata.var.columns:
            # Use genes that are variable in any batch
            adata = adata[:, adata.var.highly_variable_nbatches > 0]
        else:
            adata = adata[:, adata.var.highly_variable]
    
    # Scale data if requested
    if scale_data:
        logger.info("Scaling data")
        sc.pp.scale(adata)
    
    # Calculate processing time
    end_time = datetime.now()
    processing_time = (end_time - start_time).total_seconds()
    
    # Create results
    results = {
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'n_hvg': adata.var.highly_variable.sum() if 'highly_variable' in adata.var.columns else 0,
        'processing_time_seconds': processing_time
    }
    
    return adata, results

def main():
    """Main entry point for the normalization batch job"""
    try:
        # Update job status to running
        update_job_status('RUNNING')
        
        # Get input data
        adata = get_input_data()
        
        # Get job parameters
        params = get_job_params()
        
        # Process data
        processed_adata, results = process_normalization(adata, params)
        
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