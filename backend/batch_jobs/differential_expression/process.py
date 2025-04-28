#!/usr/bin/env python
"""
Batch job for differential expression analysis
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
logger = logging.getLogger('differential_expression')

def process_differential_expression(adata, params):
    """
    Process differential expression analysis
    
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
        Results including DE genes
    """
    logger.info(f"Starting differential expression with parameters: {params}")
    
    # Get parameters
    groupby = params.get('groupby', 'clusters')
    method = params.get('method', 'wilcoxon')
    groups = params.get('groups', 'all')
    reference = params.get('reference', 'rest')
    
    # Validate parameters
    if groupby not in adata.obs.columns:
        logger.error(f"Groupby column {groupby} not found in data")
        raise ValueError(f"Groupby column {groupby} not found in data")
    
    # Start timing
    start_time = datetime.now()
    
    # If groups is 'all', use all available groups
    if groups == 'all':
        groups = adata.obs[groupby].unique().tolist()
    
    logger.info(f"Running differential expression for {len(groups)} groups")
    
    # Run differential expression
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=groups,
        reference=reference,
        method=method,
        pts=True  # Include percentage of cells expressing genes
    )
    
    # Extract results to a more portable format
    de_results = {}
    
    # Get the results
    result_keys = ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges', 'pts', 'pts_rest']
    
    for group in groups:
        de_results[group] = {}
        for key in result_keys:
            if key in adata.uns['rank_genes_groups']:
                try:
                    # Convert to Python native types for JSON serialization
                    values = adata.uns['rank_genes_groups'][key][group]
                    if hasattr(values, 'tolist'):
                        de_results[group][key] = values.tolist()
                    else:
                        de_results[group][key] = list(values)
                except Exception as e:
                    logger.warning(f"Could not extract {key} for {group}: {e}")
    
    # Calculate processing time
    end_time = datetime.now()
    processing_time = (end_time - start_time).total_seconds()
    
    # Create results
    results = {
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'n_groups': len(groups),
        'groupby': groupby,
        'method': method,
        'reference': reference,
        'processing_time_seconds': processing_time,
        'de_results': de_results
    }
    
    return adata, results

def main():
    """Main entry point for the differential expression batch job"""
    try:
        # Update job status to running
        update_job_status('RUNNING')
        
        # Get input data
        adata = get_input_data()
        
        # Get job parameters
        params = get_job_params()
        
        # Process data
        processed_adata, results = process_differential_expression(adata, params)
        
        # Save output data
        output_uris = save_output_data(processed_adata, {'de_results': results['de_results']})
        
        # Update job status to succeeded
        update_job_status('SUCCEEDED', {
            'output_s3_uri': output_uris['output_s3_uri'],
            'de_results_s3_uri': output_uris.get('de_results_s3_uri')
        })
        
        logger.info("Job completed successfully")
        return 0
    
    except Exception as e:
        logger.error(f"Job failed: {e}", exc_info=True)
        update_job_status('FAILED', {'error': str(e)})
        return 1

if __name__ == "__main__":
    sys.exit(main())