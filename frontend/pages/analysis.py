import streamlit as st
import pandas as pd
import numpy as np
import time
import json
import io

from utils.aws_utils import download_from_s3, upload_to_s3, RESULT_PREFIX
from utils.ui_components import render_data_info, render_job_status, show_job_result_summary
from utils.job_management import submit_job, check_job_status, get_job_results, process_local_job

# Try to import scanpy for data processing
try:
    import scanpy as sc
    import matplotlib.pyplot as plt
    import seaborn as sns
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False
    print("Warning: scanpy not available. Using mock data.")

def render_analysis_page():
    """Render the analysis page with job queue integration"""
    st.header("🧪 Analysis & Results")
    
    # Check if we have file info and data
    if 'file_info' not in st.session_state or not st.session_state.file_info.get('original_name'):
        st.warning("No data uploaded. Please go back to upload step.")
        if st.button("Back to Upload"):
            st.session_state.step = 0
            st.rerun()
        return
    
    # If we already have data in the session state, use it
    if st.session_state.data is not None:
        adata = st.session_state.data
    # Otherwise, try to load from S3
    elif 'filtered_s3_uri' in st.session_state.file_info:
        with st.spinner("Loading filtered data from S3..."):
            if SCANPY_AVAILABLE:
                adata = download_from_s3(st.session_state.file_info['filtered_s3_uri'])
                if adata is None:
                    st.error("Failed to load filtered data from S3.")
                    return
                st.session_state.data = adata
            else:
                # Create mock data for demonstration
                time.sleep(1)
                st.info("Using mock data since scanpy is not available.")
                adata = type('MockAnnData', (), {
                    'n_obs': 850,  # Reduced after filtering
                    'n_vars': 19000,  # Reduced after filtering
                    'obs': pd.DataFrame(index=range(850)),
                    'var': pd.DataFrame(index=[f'gene_{i}' for i in range(19000)]),
                    'var_names': [f'gene_{i}' for i in range(19000)],
                    'uns': {},
                    'obsm': {},
                    'layers': {}
                })
                st.session_state.data = adata
    # Try other data sources if filtered not available
    elif 'visualized_s3_uri' in st.session_state.file_info:
        with st.spinner("Loading visualized data from S3..."):
            adata = download_from_s3(st.session_state.file_info['visualized_s3_uri'])
            if adata is None:
                st.error("Failed to load visualized data from S3.")
                return
            st.session_state.data = adata
    elif 'prepared_s3_uri' in st.session_state.file_info:
        with st.spinner("Loading prepared data from S3..."):
            adata = download_from_s3(st.session_state.file_info['prepared_s3_uri'])
            if adata is None:
                st.error("Failed to load prepared data from S3.")
                return
            st.session_state.data = adata
    elif 's3_uri' in st.session_state.file_info:
        with st.spinner("Loading original data from S3..."):
            adata = download_from_s3(st.session_state.file_info['s3_uri'])
            if adata is None:
                st.error("Failed to load data from S3.")
                return
            st.session_state.data = adata
    else:
        st.error("No data source available. Please go back and prepare your data.")
        return
    
    # Display data info
    render_data_info(adata)
    
    st.markdown("---")
    st.markdown("### Analysis Options")
    
    # Analysis type selection
    analysis_type = st.selectbox(
        "Analysis Type",
        ["Normalization & HVG Selection", "PCA & UMAP", "Clustering", "Differential Expression"]
    )
    
    # Track if we need to show job submission/status UI
    show_job_ui = False
    job_params = {}
    
    # Prepare input data for job
    input_data = {
        'filename': st.session_state.file_info['unique_name'],
        'data_size_mb': getattr(adata, 'n_obs', 1000) * getattr(adata, 'n_vars', 20000) * 8 / 1e6  # Rough estimate of size in MB
    }
    
    # Different UI based on analysis type
    if analysis_type == "Normalization & HVG Selection":
        show_job_ui, job_params = render_normalization_ui(adata)
    elif analysis_type == "PCA & UMAP":
        show_job_ui, job_params = render_dim_reduction_ui(adata)
    elif analysis_type == "Clustering":
        show_job_ui, job_params = render_clustering_ui(adata)
    elif analysis_type == "Differential Expression":
        show_job_ui, job_params = render_differential_expression_ui(adata)
    
    # Job submission and monitoring
    if show_job_ui:
        st.markdown("---")
        st.markdown("### Job Submission")
        
        # Only show the submit button if we don't already have an active job
        if 'current_job_id' not in st.session_state or not st.session_state.current_job_id:
            if st.button("Run Analysis"):
                with st.spinner("Preparing job..."):
                    # Prepare input data
                    success, s3_uri = upload_to_s3(
                        adata,
                        f"{st.session_state.file_info['unique_name']}_analysis_input.h5ad",
                        RESULT_PREFIX
                    )
                    
                    if not success:
                        st.error("Failed to upload data to S3")
                        return
                    
                    input_data['s3_uri'] = s3_uri
                    
                    # Submit job
                    job_type = analysis_type.lower().replace(" & ", "_").replace(" ", "_")
                    job_result = submit_job(job_type, input_data, job_params)
                    
                    if 'job_id' in job_result:
                        st.session_state.current_job_id = job_result['job_id']
                        st.session_state.job_type = job_type
                        st.session_state.job_submission_time = time.time()
                        st.success(f"Job submitted successfully!")
                        
                        # If it's a local job, process it immediately
                        if job_result.get('processing_type') == 'local':
                            with st.spinner("Processing locally..."):
                                process_local_job(job_result['job_id'])
                        
                        # Rerun to show job status
                        st.rerun()
                    else:
                        st.error(f"Failed to submit job: {job_result}")
        
        # If we have an active job, show status
        if 'current_job_id' in st.session_state and st.session_state.current_job_id:
            # Check and display job status
            job_status = render_job_status(
                st.session_state.current_job_id, 
                st.session_state.job_type if 'job_type' in st.session_state else None
            )
            
            # Show result for completed jobs
            if job_status and job_status.get('status') == 'SUCCEEDED':
                st.success("Job completed successfully!")
                
                # Get results from job
                job_results = get_job_results(st.session_state.current_job_id)
                
                if job_results and 'output_s3_uri' in job_results:
                    # Store the URI for future reference
                    st.session_state.file_info['analysis_result_uri'] = job_results['output_s3_uri']
                    
                    # Load result data
                    with st.spinner("Loading results..."):
                        if SCANPY_AVAILABLE:
                            result_adata = download_from_s3(job_results['output_s3_uri'])
                            
                            if result_adata is not None:
                                # Update session state with new data
                                st.session_state.data = result_adata
                                
                                # Show result summary
                                show_job_result_summary(job_results, st.session_state.job_type)
                                
                                # Additional visualizations based on job type
                                if st.session_state.job_type == 'normalization_hvg_selection':
                                    if 'highly_variable' in result_adata.var.columns:
                                        n_hvg = result_adata.var.highly_variable.sum()
                                        
                                        # Show HVG dispersion plot
                                        try:
                                            fig = plt.figure(figsize=(10, 6))
                                            sc.pl.highly_variable_genes(result_adata, show=False)
                                            st.pyplot(fig)
                                        except Exception as e:
                                            st.error(f"Error creating HVG plot: {e}")
                                
                                elif st.session_state.job_type == 'pca_umap':
                                    # Show PCA and UMAP plots
                                    col1, col2 = st.columns(2)
                                    
                                    with col1:
                                        if 'X_pca' in result_adata.obsm:
                                            try:
                                                fig = plt.figure(figsize=(8, 8))
                                                sc.pl.pca(result_adata, show=False)
                                                st.pyplot(fig)
                                            except Exception as e:
                                                st.error(f"Error creating PCA plot: {e}")
                                    
                                    with col2:
                                        if 'X_umap' in result_adata.obsm:
                                            try:
                                                fig = plt.figure(figsize=(8, 8))
                                                sc.pl.umap(result_adata, show=False)
                                                st.pyplot(fig)
                                            except Exception as e:
                                                st.error(f"Error creating UMAP plot: {e}")
                                
                                elif st.session_state.job_type == 'clustering':
                                    # Show cluster plots
                                    if 'clusters' in result_adata.obs.columns:
                                        try:
                                            fig = plt.figure(figsize=(10, 8))
                                            if 'X_umap' in result_adata.obsm:
                                                sc.pl.umap(result_adata, color='clusters', show=False)
                                            elif 'X_pca' in result_adata.obsm:
                                                sc.pl.pca(result_adata, color='clusters', show=False)
                                            st.pyplot(fig)
                                            
                                            # Also show spatial plot if available
                                            if 'spatial' in result_adata.obsm:
                                                fig = plt.figure(figsize=(10, 10))
                                                try:
                                                    sc.pl.spatial(result_adata, color='clusters', show=False)
                                                    st.pyplot(fig)
                                                except Exception as e:
                                                    st.error(f"Error creating spatial plot: {e}")
                                        except Exception as e:
                                            st.error(f"Error creating cluster plots: {e}")
                                
                                elif st.session_state.job_type == 'differential_expression':
                                    try:
                                        # Show DE results
                                        if 'rank_genes_groups' in result_adata.uns:
                                            st.markdown("### Differential Expression Results")
                                            
                                            # If there are DE results in a separate file
                                            if 'de_results_s3_uri' in job_results:
                                                de_results = download_from_s3(job_results['de_results_s3_uri'])
                                                
                                                if de_results:
                                                    # Let user select group to view
                                                    groups = list(de_results.keys())
                                                    selected_group = st.selectbox("Select group", groups)
                                                    
                                                    if selected_group:
                                                        group_data = de_results[selected_group]
                                                        
                                                        # Create dataframe of results
                                                        if 'names' in group_data and len(group_data['names']) > 0:
                                                            df = pd.DataFrame({
                                                                'gene': group_data.get('names', [])[:20],
                                                                'log_fold_change': group_data.get('logfoldchanges', [])[:20],
                                                                'p_value': group_data.get('pvals', [])[:20],
                                                                'p_value_adj': group_data.get('pvals_adj', [])[:20]
                                                            })
                                                            
                                                            st.write(f"Top differentially expressed genes for {selected_group}:")
                                                            st.dataframe(df)
                                                            
                                                            # Plot top genes
                                                            try:
                                                                fig = plt.figure(figsize=(12, 5))
                                                                sc.pl.rank_genes_groups_violin(
                                                                    result_adata, 
                                                                    groups=[selected_group],
                                                                    n_genes=5,
                                                                    show=False
                                                                )
                                                                st.pyplot(fig)
                                                            except Exception as e:
                                                                st.error(f"Error creating gene expression plot: {e}")
                                    except Exception as e:
                                        st.error(f"Error loading differential expression results: {e}")
                            else:
                                st.error("Failed to load results data from S3")
                        else:
                            # Mock results for demonstration
                            time.sleep(1)
                            st.info("Using mock results since scanpy is not available.")
                            
                            # Show mock result based on job type
                            mock_job_results = {
                                'status': 'SUCCESS',
                                'n_cells': 850,
                                'n_genes': 19000
                            }
                            
                            if st.session_state.job_type == 'normalization_hvg_selection':
                                mock_job_results['n_hvg'] = 2000
                            elif st.session_state.job_type == 'clustering':
                                mock_job_results['n_clusters'] = 8
                            
                            show_job_result_summary(mock_job_results, st.session_state.job_type)
                else:
                    st.error("No results found for job")
                
                # Allow user to clear job and run another analysis
                if st.button("Run Another Analysis"):
                    st.session_state.current_job_id = None
                    st.session_state.job_type = None
                    st.rerun()
    
    # Navigation button
    if st.button("⬅️ Back to Quality Control"):
        st.session_state.step = 3
        st.rerun()

def render_normalization_ui(adata):
    """
    Render UI for normalization and return job parameters
    
    Returns:
    --------
    bool
        Whether to show job UI
    dict
        Job parameters
    """
    st.subheader("Normalization & HVG Selection")
    
    # Check if data is already normalized
    is_normalized = False
    if hasattr(adata, 'uns') and 'log1p' in adata.uns:
        is_normalized = True
        st.info("Data appears to be already normalized.")
    
    # Normalization options
    col1, col2 = st.columns(2)
    
    with col1:
        normalize = st.checkbox("Normalize total counts per cell", value=not is_normalized)
        log_transform = st.checkbox("Log transform", value=not is_normalized)
    
    with col2:
        n_top_genes = st.slider("Number of top variable genes", 500, 5000, 2000)
    
    # Batch processing if batch key exists
    use_batch = False
    batch_key = None
    
    if hasattr(adata, 'obs'):
        potential_batch_keys = [
            col for col in adata.obs.columns 
            if (hasattr(adata.obs[col], 'dtype') and adata.obs[col].dtype.name == 'category') or 
               len(adata.obs[col].unique()) <= 20
        ]
        
        if potential_batch_keys:
            use_batch = st.checkbox("Account for batch effects", value=False)
            
            if use_batch:
                batch_key = st.selectbox(
                    "Batch column", 
                    potential_batch_keys,
                    index=0 if 'batch' in potential_batch_keys else 0
                )
    
    # Additional options
    col1, col2 = st.columns(2)
    
    with col1:
        subset_hvg = st.checkbox("Subset to highly variable genes", value=True)
    
    with col2:
        scale_data = st.checkbox("Scale data", value=True)
    
    # Create job parameters
    job_params = {
        'normalize': normalize,
        'log_transform': log_transform,
        'n_top_genes': n_top_genes,
        'use_batch': use_batch,
        'batch_key': batch_key,
        'subset_hvg': subset_hvg,
        'scale_data': scale_data,
        'step': 'preparation'
    }
    
    return True, job_params

def render_dim_reduction_ui(adata):
    """
    Render UI for dimensionality reduction and return job parameters
    
    Returns:
    --------
    bool
        Whether to show job UI
    dict
        Job parameters
    """
    st.subheader("PCA & UMAP")
    
    # Check if HVGs have been selected
    if hasattr(adata, 'var') and 'highly_variable' in adata.var.columns:
        has_hvg = True
        n_hvgs = adata.var.highly_variable.sum() if hasattr(adata.var.highly_variable, 'sum') else 0
        st.info(f"Dataset has {n_hvgs} highly variable genes.")
    else:
        has_hvg = False
        st.warning("Highly variable genes have not been selected. Consider running normalization first.")
    
    # Options
    col1, col2 = st.columns(2)
    
    with col1:
        n_pcs = st.slider("Number of principal components", 10, 100, 50)
    
    with col2:
        n_neighbors = st.slider("Number of neighbors for UMAP", 5, 50, 15)
    
    # Use HVG option if available
    use_hvg = False
    if has_hvg:
        use_hvg = st.checkbox("Use only highly variable genes", value=True)
    
    # Processing steps
    col1, col2, col3 = st.columns(3)
    
    with col1:
        run_pca = st.checkbox("Run PCA", value=True)
    
    with col2:
        run_neighbors = st.checkbox("Compute neighborhood graph", value=True)
    
    with col3:
        run_umap = st.checkbox("Run UMAP", value=True)
    
    # Create job parameters
    job_params = {
        'n_pcs': n_pcs,
        'n_neighbors': n_neighbors,
        'use_hvg': use_hvg if has_hvg else False,
        'run_pca': run_pca,
        'run_neighbors': run_neighbors,
        'run_umap': run_umap,
        'step': 'analysis'
    }
    
    return True, job_params

def render_clustering_ui(adata):
    """
    Render UI for clustering and return job parameters
    
    Returns:
    --------
    bool
        Whether to show job UI
    dict
        Job parameters
    """
    st.subheader("Clustering")
    
    # Check if neighbors have been computed
    has_neighbors = False
    if hasattr(adata, 'uns') and 'neighbors' in adata.uns:
        has_neighbors = True
    
    if not has_neighbors:
        st.warning("Neighborhood graph has not been computed. Consider running dimensionality reduction first.")
        show_job_ui = False
    else:
        # Options
        col1, col2 = st.columns(2)
        
        with col1:
            resolution = st.slider("Clustering resolution", 0.1, 2.0, 0.5, 0.1)
        
        with col2:
            cluster_method = st.selectbox(
                "Clustering algorithm",
                ["leiden", "louvain"],
                index=0
            )
        
        key_added = st.text_input("Column name for clusters", "clusters")
        
        # Create job parameters
        job_params = {
            'resolution': resolution,
            'cluster_method': cluster_method,
            'key_added': key_added,
            'step': 'analysis'
        }
        
        return True, job_params
    
    return False, {}

def render_differential_expression_ui(adata):
    """
    Render UI for differential expression and return job parameters
    
    Returns:
    --------
    bool
        Whether to show job UI
    dict
        Job parameters
    """
    st.subheader("Differential Expression Analysis")
    
    # Check if clustering has been done
    has_clusters = False
    potential_group_cols = []
    
    if hasattr(adata, 'obs'):
        potential_group_cols = [
            col for col in adata.obs.columns 
            if (hasattr(adata.obs[col], 'dtype') and adata.obs[col].dtype.name == 'category') or 
               len(adata.obs[col].unique()) > 1 and len(adata.obs[col].unique()) < 50
        ]
        
        has_clusters = 'clusters' in adata.obs.columns or len(potential_group_cols) > 0
    
    if not has_clusters:
        st.warning("No clustering results found. Please run clustering first.")
        return False, {}
    
    # Options
    col1, col2 = st.columns(2)
    
    with col1:
        groupby = st.selectbox(
            "Group by", 
            potential_group_cols,
            index=potential_group_cols.index('clusters') if 'clusters' in potential_group_cols else 0
        )
    
    # Get groups for selected column
    if hasattr(adata, 'obs') and groupby in adata.obs.columns:
        groups = list(adata.obs[groupby].unique())
    else:
        groups = []
    
    # Let user select specific groups or all
    with col2:
        group_selection = st.radio(
            "Select groups to analyze",
            ["All groups", "Select specific groups"],
            index=0
        )
    
    selected_groups = groups
    if group_selection == "Select specific groups":
        selected_groups = st.multiselect(
            "Select groups", 
            groups,
            default=groups[:1] if groups else []
        )
        
        if not selected_groups:
            st.warning("Please select at least one group")
            return False, {}
    
    # Additional options
    col1, col2 = st.columns(2)
    
    with col1:
        method = st.selectbox(
            "Statistical method",
            ["wilcoxon", "t-test", "logreg"],
            index=0
        )
    
    with col2:
        reference = st.selectbox(
            "Reference",
            ["rest"] + groups,
            index=0
        )
    
    # Create job parameters
    job_params = {
        'groupby': groupby,
        'groups': 'all' if group_selection == "All groups" else selected_groups,
        'method': method,
        'reference': reference,
        'step': 'analysis'
    }
    
    return True, job_params

if __name__ == "__main__":
    # For testing this page in isolation
    render_analysis_page()