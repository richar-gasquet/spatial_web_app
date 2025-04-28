import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time

def render_header():
    """Render application header with title and description"""
    st.title("📊 Spatial Transcriptomics Analysis Pipeline")
    
    # Only show description on the first page
    if st.session_state.step == 0:
        st.markdown("""
        This application provides a multi-step workflow for analyzing spatial transcriptomics data:
        
        1. **Upload Data** - Upload your spatial data files
        2. **Data Preparation** - Calculate QC metrics and prepare data
        3. **Data Visualization** - Visualize data distribution and QC metrics
        4. **Quality Control & Filtering** - Filter spots/cells and genes
        5. **Analysis** - Run various analyses including normalization, clustering, and differential expression
        
        Large or computationally intensive tasks will automatically be processed in the cloud.
        """)

def render_sidebar():
    """Render sidebar navigation"""
    st.sidebar.title("Navigation")
    
    steps = [
        "1. Upload Data",
        "2. Data Preparation",
        "3. Data Visualization",
        "4. Quality Control & Filtering",
        "5. Analysis & Results"
    ]
    
    # Show navigation (but only allow clicking on steps that are available)
    for i, step_name in enumerate(steps):
        if i <= st.session_state.step:  # Can only navigate to current or previous steps
            if st.sidebar.button(step_name, key=f"nav_{i}"):
                st.session_state.step = i
                st.rerun()
        else:
            # Show disabled button
            st.sidebar.button(step_name, disabled=True, key=f"nav_{i}")
    
    # Display progress
    st.sidebar.progress((st.session_state.step) / (len(steps) - 1))
    
    # Start over button
    if st.sidebar.button("↩️ Start Over"):
        for key in list(st.session_state.keys()):
            if key not in ['aws_clients']:  # Preserve AWS clients
                del st.session_state[key]
        st.session_state.step = 0
        st.rerun()
    
    # Display session info
    st.sidebar.markdown("---")
    if 'file_info' in st.session_state and st.session_state.file_info.get('original_name'):
        st.sidebar.markdown(f"**Current file:** {st.session_state.file_info['original_name']}")
    
    # AWS connection status
    if 'aws_clients' in st.session_state:
        if 'error' in st.session_state.aws_clients:
            st.sidebar.error("⚠️ AWS connection failed")
        else:
            st.sidebar.success("✅ AWS connected")
    else:
        st.sidebar.warning("⚠️ AWS not connected")

def render_job_status(job_id, job_type):
    """
    Render job status with automatic refresh
    
    Parameters:
    -----------
    job_id : str
        Job ID to check
    job_type : str
        Type of job (for display purposes)
    
    Returns:
    --------
    dict
        Latest job status
    """
    from utils.job_management import check_job_status
    
    # If not provided, check session state
    if not job_id and 'current_job_id' in st.session_state:
        job_id = st.session_state.current_job_id
    
    if not job_type and 'job_type' in st.session_state:
        job_type = st.session_state.job_type
    
    if not job_id:
        return None
    
    # Set initial submission time if not set
    if 'job_submission_time' not in st.session_state:
        st.session_state.job_submission_time = time.time()
    
    # Poll for job status
    job_status = check_job_status(job_id)
    
    # Display job status
    status_color = {
        'SUBMITTED': 'blue',
        'SUBMITTED_TO_BATCH': 'blue',
        'RUNNING': 'orange',
        'SUCCEEDED': 'green',
        'FAILED': 'red'
    }
    
    # Calculate elapsed time
    elapsed_time = time.time() - st.session_state.job_submission_time
    elapsed_min = int(elapsed_time // 60)
    elapsed_sec = int(elapsed_time % 60)
    
    # Show status with color
    status = job_status.get('status', 'UNKNOWN')
    
    # Create a status container
    status_container = st.empty()
    
    # Show status with appropriate styling
    status_container.markdown(
        f"**Job Status:** <span style='color:{status_color.get(status, 'gray')}'>{status}</span> "
        f"(Running for {elapsed_min}m {elapsed_sec}s)", 
        unsafe_allow_html=True
    )
    
    # Show progress bar for running jobs
    if status in ['SUBMITTED', 'SUBMITTED_TO_BATCH', 'RUNNING']:
        progress_bar = st.progress(0.5)  # Indeterminate progress
        
        # Add a refresh button
        refresh_col, auto_col = st.columns([1, 2])
        with refresh_col:
            if st.button("🔄 Refresh"):
                return check_job_status(job_id)  # Immediate refresh
        
        with auto_col:
            st.info("Status updates automatically every 30 seconds")
        
        # Schedule a rerun after a delay for auto-refresh
        if 'last_refresh' not in st.session_state:
            st.session_state.last_refresh = 0
        
        if time.time() - st.session_state.last_refresh > 30:
            st.session_state.last_refresh = time.time()
            time.sleep(0.1)  # Small delay
            st.rerun()
    
    # Show error for failed jobs
    elif status == 'FAILED':
        st.error("Job failed")
        
        # Show error message if available
        if 'error' in job_status:
            st.error(f"Error: {job_status['error']}")
    
    return job_status

def render_data_info(adata):
    """
    Render basic information about the dataset
    
    Parameters:
    -----------
    adata : AnnData
        The dataset to show info for
    """
    if adata is None:
        st.warning("No data available")
        return
    
    # Basic dataset info
    st.markdown("### Dataset Information")
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Number of spots/cells", adata.n_obs)
    with col2:
        st.metric("Number of genes", adata.n_vars)
    
    # Additional metadata if available
    if hasattr(adata, 'obs') and len(adata.obs.columns) > 0:
        st.markdown("### Metadata")
        
        # Show first few rows of metadata
        metadata_sample = adata.obs.head(5)
        st.dataframe(metadata_sample)
        
        # Show categorical columns as options for grouping
        cat_cols = [
            col for col in adata.obs.columns 
            if adata.obs[col].dtype.name == 'category' or 
            len(adata.obs[col].unique()) < 20
        ]
        
        if cat_cols:
            st.markdown("**Categorical columns available:**")
            st.write(", ".join(cat_cols))

def show_job_result_summary(job_result, job_type):
    """
    Show a summary of job results based on job type
    
    Parameters:
    -----------
    job_result : dict
        Job result dictionary
    job_type : str
        Type of job
    """
    if not job_result or 'status' not in job_result or job_result['status'] != 'SUCCESS':
        st.error("No valid results to display")
        return
    
    st.markdown("### Analysis Results")
    
    # Display different summaries based on job type
    if job_type == 'normalization':
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Cells/Spots", job_result.get('n_cells', 'N/A'))
        with col2:
            st.metric("Genes", job_result.get('n_genes', 'N/A')) 
        with col3:
            st.metric("Highly Variable Genes", job_result.get('n_hvg', 'N/A'))
            
    elif job_type == 'dimensionality_reduction':
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Cells/Spots", job_result.get('n_cells', 'N/A'))
        with col2:
            st.metric("Genes", job_result.get('n_genes', 'N/A'))
        
        # Show dimensionality reduction info if available
        if 'pca_shape' in job_result:
            st.write(f"PCA representation shape: {job_result['pca_shape']}")
        if 'umap_shape' in job_result:
            st.write(f"UMAP representation shape: {job_result['umap_shape']}")
            
    elif job_type == 'clustering':
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Cells/Spots", job_result.get('n_cells', 'N/A'))
        with col2:
            st.metric("Genes", job_result.get('n_genes', 'N/A'))
        with col3:
            st.metric("Clusters", job_result.get('n_clusters', 'N/A'))
            
        # Show cluster distribution if available
        if 'cluster_counts' in job_result:
            cluster_df = pd.DataFrame(
                list(job_result['cluster_counts'].items()),
                columns=['Cluster', 'Count']
            )
            cluster_df = cluster_df.sort_values('Cluster')
            
            st.markdown("**Cluster distribution:**")
            st.dataframe(cluster_df)
            
    elif job_type == 'differential_expression':
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Cells/Spots", job_result.get('n_cells', 'N/A'))
        with col2:
            st.metric("Genes", job_result.get('n_genes', 'N/A'))
        with col3:
            st.metric("Groups", job_result.get('n_groups', 'N/A'))
            
        st.markdown("Results include differentially expressed genes for each group.")
        st.info("To view the full results, load the dataset and use the visualization options.")