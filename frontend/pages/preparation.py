import streamlit as st
import pandas as pd
import numpy as np
import time
import json

from utils.aws_utils import download_from_s3, upload_to_s3, PREPARED_PREFIX
from utils.ui_components import render_data_info
from utils.job_management import submit_job, check_job_status, get_job_results

# Try to import scanpy for data processing
try:
    import scanpy as sc
    import anndata as ad
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False
    print("Warning: scanpy not available. Using mock data.")

def render_preparation_page():
    """Render the data preparation page"""
    st.header("🔍 Data Preparation")
    
    # Check if we have file info
    if 'file_info' not in st.session_state or not st.session_state.file_info.get('original_name'):
        st.warning("No data uploaded. Please go back to upload step.")
        if st.button("Back to Upload"):
            st.session_state.step = 0
            st.rerun()
        return
    
    # Display file info
    st.write(f"File: **{st.session_state.file_info['original_name']}**")
    
    # If we already have data in the session state, use it
    if st.session_state.data is not None:
        adata = st.session_state.data
    # Otherwise, try to load from S3
    elif 's3_uri' in st.session_state.file_info:
        with st.spinner("Loading data from S3..."):
            if SCANPY_AVAILABLE:
                # Load actual data
                adata = download_from_s3(st.session_state.file_info['s3_uri'])
                if adata is None:
                    st.error("Failed to load data from S3.")
                    return
                st.session_state.data = adata
            else:
                # Create mock data for demonstration
                time.sleep(1)  # Simulate loading
                st.info("Using mock data since scanpy is not available.")
                # Create minimal fake AnnData object
                adata = type('MockAnnData', (), {
                    'n_obs': 1000,
                    'n_vars': 20000,
                    'obs': pd.DataFrame(index=range(1000)),
                    'var': pd.DataFrame(index=[f'gene_{i}' for i in range(20000)]),
                    'var_names': [f'gene_{i}' for i in range(20000)],
                    'uns': {},
                    'obsm': {},
                    'layers': {}
                })
                st.session_state.data = adata
    else:
        st.error("No data source available. Please re-upload your file.")
        return
    
    # Display data info
    render_data_info(adata)
    
    st.markdown("---")
    st.markdown("### Preparation Options")
    
    # QC metrics calculation
    st.subheader("Quality Control Metrics")
    
    col1, col2 = st.columns(2)
    
    with col1:
        add_mt = st.checkbox("Add mitochondrial gene annotation", value=True)
        add_hb = st.checkbox("Add hemoglobin gene annotation", value=True)
    
    with col2:
        calculate_qc = st.checkbox("Calculate QC metrics", value=True)
    
    # Process data button
    if st.button("Calculate QC Metrics"):
        with st.spinner("Calculating QC metrics..."):
            if SCANPY_AVAILABLE:
                # Perform actual QC calculation
                
                # Add gene annotations if requested
                if add_mt:
                    adata.var['mt'] = adata.var_names.str.startswith(('mt-', 'MT-'))
                    st.write("Added mitochondrial gene annotations.")
                
                if add_hb:
                    adata.var['hb'] = adata.var_names.str.contains(('Hb.*-', 'HB.*-'))
                    st.write("Added hemoglobin gene annotations.")
                
                # Calculate QC metrics
                if calculate_qc:
                    qc_vars = []
                    if add_mt:
                        qc_vars.append('mt')
                    if add_hb:
                        qc_vars.append('hb')
                    
                    sc.pp.calculate_qc_metrics(
                        adata, 
                        qc_vars=qc_vars if qc_vars else None, 
                        percent_top=None, 
                        log1p=False, 
                        inplace=True
                    )
                    st.write("Calculated QC metrics.")
                
                # Update session state
                st.session_state.data = adata
                
                # Save to S3
                with st.spinner("Saving prepared data to S3..."):
                    success, s3_uri = upload_to_s3(
                        adata,
                        f"{st.session_state.file_info['unique_name']}_prepared.h5ad",
                        PREPARED_PREFIX
                    )
                    
                    if success:
                        st.success("Data preparation complete! Saved to S3.")
                        st.session_state.file_info['prepared_s3_uri'] = s3_uri
                    else:
                        st.error("Failed to save prepared data to S3.")
            else:
                # Mock processing for demonstration
                time.sleep(2)  # Simulate processing time
                st.success("QC metrics calculated (mock data).")
                
                # Add mock QC metrics to our fake data
                adata.obs['n_genes_by_counts'] = np.random.randint(500, 5000, size=adata.n_obs)
                adata.obs['total_counts'] = adata.obs['n_genes_by_counts'] * np.random.randint(5, 50, size=adata.n_obs)
                adata.obs['pct_counts_mt'] = np.random.uniform(0, 30, size=adata.n_obs)
                
                # Update session state
                st.session_state.data = adata
    
    # Show QC metrics if available
    if hasattr(adata, 'obs') and 'n_genes_by_counts' in adata.obs.columns:
        st.markdown("### QC Metrics Preview")
        
        # Show QC metrics table
        qc_cols = [col for col in adata.obs.columns if any(x in col for x in ['counts', 'pct'])]
        if qc_cols:
            st.dataframe(adata.obs[qc_cols].head())
            
            # Show histograms
            st.markdown("### QC Metrics Distribution")
            
            if SCANPY_AVAILABLE:
                # Real plots with scanpy
                try:
                    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
                    
                    sns.histplot(adata.obs['n_genes_by_counts'], kde=True, ax=axs[0])
                    axs[0].set_title('Genes per Cell')
                    
                    sns.histplot(adata.obs['total_counts'], kde=True, ax=axs[1])
                    axs[1].set_title('UMI Counts per Cell')
                    
                    if 'pct_counts_mt' in adata.obs.columns:
                        sns.histplot(adata.obs['pct_counts_mt'], kde=True, ax=axs[2])
                        axs[2].set_title('% Mitochondrial')
                    
                    plt.tight_layout()
                    st.pyplot(fig)
                except Exception as e:
                    st.error(f"Error creating QC plots: {e}")
            else:
                # Mock plots
                st.info("QC plots would be shown here. Install matplotlib and seaborn to see actual plots.")
    
    # Navigation buttons
    col1, col2 = st.columns(2)
    
    with col1:
        if st.button("⬅️ Back to Upload"):
            st.session_state.step = 0
            st.rerun()
    
    with col2:
        if st.button("Proceed to Visualization ➡️"):
            st.session_state.step = 2
            st.rerun()

if __name__ == "__main__":
    # For testing this page in isolation
    render_preparation_page()