import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import io
import tempfile
import os

from utils.aws_utils import upload_to_s3, download_from_s3, RESULT_PREFIX

# Import scanpy only if available
try:
    import scanpy as sc
    import anndata as ad
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False
    print("Warning: scanpy not installed. Analysis functions will not be available.")

def load_data_from_s3(s3_uri):
    """Load data from S3 and return AnnData object"""
    # Use AWS utils to download
    return download_from_s3(s3_uri)

def process_normalization(input_data, job_params):
    """
    Process normalization and HVG selection
    
    Parameters:
    -----------
    input_data : dict
        Contains input S3 URI
    job_params : dict
        Parameters for normalization
        
    Returns:
    --------
    dict
        Results including output S3 URI
    """
    if not SCANPY_AVAILABLE:
        return {
            'status': 'ERROR',
            'error': 'scanpy not installed'
        }
    
    # Load data from S3
    adata = load_data_from_s3(input_data['s3_uri'])
    
    if adata is None:
        return {
            'status': 'ERROR',
            'error': 'Failed to load data from S3'
        }
    
    try:
        # Get parameters
        normalize = job_params.get('normalize', True)
        log_transform = job_params.get('log_transform', True)
        n_top_genes = job_params.get('n_top_genes', 2000)
        use_batch = job_params.get('use_batch', False)
        batch_key = job_params.get('batch_key', None)
        scale_data = job_params.get('scale_data', True)
        
        # Make a copy of counts for later use if needed
        if 'counts' not in adata.layers:
            adata.layers['counts'] = adata.X.copy()
        
        # Normalize
        if normalize:
            sc.pp.normalize_total(adata, inplace=True)
        
        # Log transform
        if log_transform:
            sc.pp.log1p(adata)
        
        # Find highly variable genes
        if use_batch and batch_key in adata.obs.columns:
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
        adata.raw = adata
        
        # Subset to highly variable genes if requested
        if job_params.get('subset_hvg', True):
            # If batch mode, use genes that are variable in any batch
            if use_batch and batch_key in adata.obs.columns:
                adata = adata[:, adata.var.highly_variable_nbatches > 0]
            else:
                adata = adata[:, adata.var.highly_variable]
        
        # Scale data if requested
        if scale_data:
            sc.pp.scale(adata)
        
        # Save processed data to S3
        output_filename = f"{input_data['filename']}_normalized.h5ad"
        success, output_s3_uri = upload_to_s3(
            adata,
            output_filename,
            RESULT_PREFIX
        )
        
        if not success:
            return {
                'status': 'ERROR',
                'error': 'Failed to upload processed data to S3'
            }
        
        # Return results
        return {
            'status': 'SUCCESS',
            'output_s3_uri': output_s3_uri,
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'n_hvg': adata.var.highly_variable.sum() if 'highly_variable' in adata.var.columns else 0
        }
    
    except Exception as e:
        return {
            'status': 'ERROR',
            'error': str(e)
        }

def process_dimensionality_reduction(input_data, job_params):
    """
    Process dimensionality reduction (PCA, neighbors, UMAP)
    
    Parameters:
    -----------
    input_data : dict
        Contains input S3 URI
    job_params : dict
        Parameters for dimensionality reduction
        
    Returns:
    --------
    dict
        Results including output S3 URI
    """
    if not SCANPY_AVAILABLE:
        return {
            'status': 'ERROR',
            'error': 'scanpy not installed'
        }
    
    # Load data from S3
    adata = load_data_from_s3(input_data['s3_uri'])
    
    if adata is None:
        return {
            'status': 'ERROR',
            'error': 'Failed to load data from S3'
        }
    
    try:
        # Get parameters
        n_pcs = job_params.get('n_pcs', 50)
        n_neighbors = job_params.get('n_neighbors', 15)
        run_pca = job_params.get('run_pca', True)
        run_neighbors = job_params.get('run_neighbors', True)
        run_umap = job_params.get('run_umap', True)
        
        # Run PCA
        if run_pca:
            sc.tl.pca(adata, n_comps=n_pcs)
        
        # Run neighbors
        if run_neighbors:
            sc.pp.neighbors(adata, n_neighbors=n_neighbors)
        
        # Run UMAP
        if run_umap and 'neighbors' in adata.uns:
            sc.tl.umap(adata)
        
        # Save processed data to S3
        output_filename = f"{input_data['filename']}_dim_reduced.h5ad"
        success, output_s3_uri = upload_to_s3(
            adata,
            output_filename,
            RESULT_PREFIX
        )
        
        if not success:
            return {
                'status': 'ERROR',
                'error': 'Failed to upload processed data to S3'
            }
        
        # Return results
        result = {
            'status': 'SUCCESS',
            'output_s3_uri': output_s3_uri,
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars
        }
        
        # Add info about computed representations
        if 'X_pca' in adata.obsm:
            result['pca_shape'] = adata.obsm['X_pca'].shape
        
        if 'X_umap' in adata.obsm:
            result['umap_shape'] = adata.obsm['X_umap'].shape
        
        return result
    
    except Exception as e:
        return {
            'status': 'ERROR',
            'error': str(e)
        }

def process_clustering(input_data, job_params):
    """
    Process clustering
    
    Parameters:
    -----------
    input_data : dict
        Contains input S3 URI
    job_params : dict
        Parameters for clustering
        
    Returns:
    --------
    dict
        Results including output S3 URI
    """
    if not SCANPY_AVAILABLE:
        return {
            'status': 'ERROR',
            'error': 'scanpy not installed'
        }
    
    # Load data from S3
    adata = load_data_from_s3(input_data['s3_uri'])
    
    if adata is None:
        return {
            'status': 'ERROR',
            'error': 'Failed to load data from S3'
        }
    
    try:
        # Get parameters
        resolution = job_params.get('resolution', 0.5)
        cluster_method = job_params.get('cluster_method', 'leiden')
        key_added = job_params.get('key_added', 'clusters')
        
        # Check if neighbors have been computed
        if 'neighbors' not in adata.uns:
            return {
                'status': 'ERROR',
                'error': 'Neighborhood graph not found. Run dimensionality reduction first.'
            }
        
        # Run clustering
        if cluster_method.lower() == 'leiden':
            sc.tl.leiden(adata, resolution=resolution, key_added=key_added)
        elif cluster_method.lower() == 'louvain':
            sc.tl.louvain(adata, resolution=resolution, key_added=key_added)
        else:
            return {
                'status': 'ERROR',
                'error': f'Unknown clustering method: {cluster_method}'
            }
        
        # Save processed data to S3
        output_filename = f"{input_data['filename']}_clustered.h5ad"
        success, output_s3_uri = upload_to_s3(
            adata,
            output_filename,
            RESULT_PREFIX
        )
        
        if not success:
            return {
                'status': 'ERROR',
                'error': 'Failed to upload processed data to S3'
            }
        
        # Return results
        return {
            'status': 'SUCCESS',
            'output_s3_uri': output_s3_uri,
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'n_clusters': len(adata.obs[key_added].unique()),
            'cluster_counts': adata.obs[key_added].value_counts().to_dict()
        }
    
    except Exception as e:
        return {
            'status': 'ERROR',
            'error': str(e)
        }

def process_differential_expression(input_data, job_params):
    """
    Process differential expression analysis
    
    Parameters:
    -----------
    input_data : dict
        Contains input S3 URI
    job_params : dict
        Parameters for differential expression
        
    Returns:
    --------
    dict
        Results including output S3 URI
    """
    if not SCANPY_AVAILABLE:
        return {
            'status': 'ERROR',
            'error': 'scanpy not installed'
        }
    
    # Load data from S3
    adata = load_data_from_s3(input_data['s3_uri'])
    
    if adata is None:
        return {
            'status': 'ERROR',
            'error': 'Failed to load data from S3'
        }
    
    try:
        # Get parameters
        groupby = job_params.get('groupby', 'clusters')
        method = job_params.get('method', 'wilcoxon')
        groups = job_params.get('groups', 'all')
        reference = job_params.get('reference', 'rest')
        
        # Validate parameters
        if groupby not in adata.obs.columns:
            return {
                'status': 'ERROR',
                'error': f'Groupby column {groupby} not found in data'
            }
        
        # If groups is 'all', use all available groups
        if groups == 'all':
            groups = adata.obs[groupby].unique().tolist()
        
        # Run differential expression
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            groups=groups,
            reference=reference,
            method=method
        )
        
        # Extract results to a more portable format
        de_results = {}
        
        # Get the results
        result_keys = ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']
        for group in groups:
            de_results[group] = {}
            for key in result_keys:
                if key in adata.uns['rank_genes_groups']:
                    try:
                        de_results[group][key] = adata.uns['rank_genes_groups'][key][group].tolist()
                    except:
                        # If not convertible to list, skip
                        pass
        
        # Save processed data to S3
        output_filename = f"{input_data['filename']}_de.h5ad"
        success, output_s3_uri = upload_to_s3(
            adata,
            output_filename,
            RESULT_PREFIX
        )
        
        if not success:
            return {
                'status': 'ERROR',
                'error': 'Failed to upload processed data to S3'
            }
        
        # Also save DE results as a separate portable file
        de_results_filename = f"{input_data['filename']}_de_results.json"
        success, de_results_s3_uri = upload_to_s3(
            de_results,
            de_results_filename,
            RESULT_PREFIX
        )
        
        # Return results
        return {
            'status': 'SUCCESS',
            'output_s3_uri': output_s3_uri,
            'de_results_s3_uri': de_results_s3_uri,
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'n_groups': len(groups)
        }
    
    except Exception as e:
        return {
            'status': 'ERROR',
            'error': str(e)
        }