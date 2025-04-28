import streamlit as st
import pandas as pd
import numpy as np
import time

from utils.aws_utils import download_from_s3, upload_to_s3, QC_PREFIX
from utils.ui_components import render_data_info

# Try to import scanpy for data processing
try:
    import scanpy as sc
    import matplotlib.pyplot as plt
    import seaborn as sns
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False
    print("Warning: scanpy not available. Using mock data.")

def render_qc_filtering_page():
    """Render the quality control and filtering page"""
    st.header("🧹 Quality Control & Filtering")
    
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
    
    # Display current data info
    render_data_info(adata)
    
    st.markdown("---")
    st.markdown("### Quality Control Filtering")
    
    # Check if QC metrics are available
    has_qc_metrics = False
    if hasattr(adata, 'obs'):
        qc_cols = [col for col in adata.obs.columns if any(x in col for x in ['counts', 'pct', 'n_genes'])]
        has_qc_metrics = len(qc_cols) > 0
    
    if not has_qc_metrics:
        st.warning("No QC metrics found. Please go back to the preparation step to calculate QC metrics.")
    else:
        # Cell filtering
        st.subheader("Cell/Spot Filtering")
        
        col1, col2 = st.columns(2)
        
        # Get reasonable defaults based on data
        if 'n_genes_by_counts' in adata.obs.columns:
            min_genes_default = int(np.percentile(adata.obs['n_genes_by_counts'], 5))
        else:
            min_genes_default = 200
            
        if 'pct_counts_mt' in adata.obs.columns:
            max_mt_default = min(30.0, np.percentile(adata.obs['pct_counts_mt'], 95))
        else:
            max_mt_default = 20.0
            
        if 'pct_counts_hb' in adata.obs.columns:
            max_hb_default = min(10.0, np.percentile(adata.obs['pct_counts_hb'], 95))
        else:
            max_hb_default = 5.0
        
        with col1:
            min_genes = st.slider(
                "Minimum genes per cell/spot", 
                0, 5000, min_genes_default,
                help="Filter out cells with fewer than this many genes detected"
            )
            
            max_mt_pct = st.slider(
                "Maximum mitochondrial %", 
                0.0, 100.0, max_mt_default,
                help="Filter out cells with higher mitochondrial content"
            )
        
        with col2:
            min_counts = st.slider(
                "Minimum counts per cell/spot", 
                0, 10000, min_genes_default * 2,
                help="Filter out cells with fewer than this many counts"
            )
            
            max_hb_pct = st.slider(
                "Maximum hemoglobin %", 
                0.0, 100.0, max_hb_default,
                help="Filter out cells with higher hemoglobin content"
            )
        
        # Gene filtering
        st.subheader("Gene Filtering")
        
        col1, col2 = st.columns(2)
        
        with col1:
            filter_mt = st.checkbox("Filter mitochondrial genes", value=True)
            filter_hb = st.checkbox("Filter hemoglobin genes", value=True)
        
        with col2:
            min_cells = st.slider(
                "Minimum cells per gene", 
                0, 100, 5,
                help="Filter out genes detected in fewer than this percentage of cells"
            )
            
            filter_custom = st.checkbox("Filter custom genes", value=False)
            
            if filter_custom:
                custom_genes = st.text_input(
                    "Enter gene names to filter (comma-separated)",
                    help="Specify genes to exclude by name"
                )
        
        # Preview filter effects
        if st.button("Preview Filter Effects"):
            with st.spinner("Calculating filter effects..."):
                if SCANPY_AVAILABLE:
                    # Construct filter expression
                    cell_filters = []
                    
                    if 'n_genes_by_counts' in adata.obs.columns:
                        cell_filters.append(f"adata.obs['n_genes_by_counts'] >= {min_genes}")
                    
                    if 'total_counts' in adata.obs.columns:
                        cell_filters.append(f"adata.obs['total_counts'] >= {min_counts}")
                    
                    if 'pct_counts_mt' in adata.obs.columns:
                        cell_filters.append(f"adata.obs['pct_counts_mt'] <= {max_mt_pct}")
                    
                    if 'pct_counts_hb' in adata.obs.columns:
                        cell_filters.append(f"adata.obs['pct_counts_hb'] <= {max_hb_pct}")
                    
                    # Count cells that would be kept
                    if cell_filters:
                        cell_filter_expr = " & ".join(cell_filters)
                        keep_cells = eval(cell_filter_expr)
                        n_cells_kept = keep_cells.sum()
                        n_cells_removed = len(keep_cells) - n_cells_kept
                    else:
                        n_cells_kept = adata.n_obs
                        n_cells_removed = 0
                    
                    # Count genes that would be kept
                    gene_filters = []
                    
                    if filter_mt and 'mt' in adata.var.columns:
                        gene_filters.append("~adata.var['mt']")
                    
                    if filter_hb and 'hb' in adata.var.columns:
                        gene_filters.append("~adata.var['hb']")
                    
                    # Min cells per gene
                    if min_cells > 0:
                        # Calculate number of cells expressing each gene
                        n_cells_per_gene = np.sum(adata.X > 0, axis=0)
                        gene_expressed_in_min_cells = n_cells_per_gene >= (min_cells/100.0 * adata.n_obs)
                        adata.var['n_cells'] = n_cells_per_gene
                        gene_filters.append(f"adata.var['n_cells'] >= {min_cells/100.0 * adata.n_obs}")
                    
                    # Custom genes to filter
                    if filter_custom and custom_genes:
                        gene_list = [g.strip() for g in custom_genes.split(',')]
                        gene_filters.append(f"~adata.var_names.isin({gene_list})")
                    
                    if gene_filters:
                        gene_filter_expr = " & ".join(gene_filters)
                        keep_genes = eval(gene_filter_expr)
                        n_genes_kept = keep_genes.sum()
                        n_genes_removed = len(keep_genes) - n_genes_kept
                    else:
                        n_genes_kept = adata.n_vars
                        n_genes_removed = 0
                    
                    # Display preview
                    st.write(f"Cells/spots that would be removed: {n_cells_removed} ({n_cells_removed/adata.n_obs:.1%})")
                    st.write(f"Genes that would be removed: {n_genes_removed} ({n_genes_removed/adata.n_vars:.1%})")
                    
                    # Show distribution plots if matplotlib is available
                    try:
                        # Create a figure with subplots
                        fig, axs = plt.subplots(1, 2, figsize=(12, 5))
                        
                        # Cell filter distribution
                        if 'n_genes_by_counts' in adata.obs.columns:
                            sns.histplot(adata.obs['n_genes_by_counts'], ax=axs[0], kde=True)
                            axs[0].axvline(x=min_genes, color='r', linestyle='--')
                            axs[0].set_title('Genes per Cell Distribution')
                            axs[0].set_xlabel('Genes per Cell')
                            axs[0].set_ylabel('Number of Cells')
                        
                        # MT content distribution
                        if 'pct_counts_mt' in adata.obs.columns:
                            sns.histplot(adata.obs['pct_counts_mt'], ax=axs[1], kde=True)
                            axs[1].axvline(x=max_mt_pct, color='r', linestyle='--')
                            axs[1].set_title('Mitochondrial Content Distribution')
                            axs[1].set_xlabel('MT Content (%)')
                            axs[1].set_ylabel('Number of Cells')
                        
                        plt.tight_layout()
                        st.pyplot(fig)
                    except Exception as e:
                        st.error(f"Error creating preview plots: {e}")
                else:
                    # Mock filter preview for demonstration
                    time.sleep(1)
                    
                    # Make up some reasonable numbers
                    total_cells = adata.n_obs
                    total_genes = adata.n_vars
                    
                    cells_removed = int(total_cells * 0.15)  # 15% cells removed
                    genes_removed = int(total_genes * 0.05)  # 5% genes removed
                    
                    st.write(f"Cells/spots that would be removed: {cells_removed} ({cells_removed/total_cells:.1%})")
                    st.write(f"Genes that would be removed: {genes_removed} ({genes_removed/total_genes:.1%})")
                    
                    st.info("Distribution plots would be shown here. Install matplotlib and seaborn to see actual plots.")
        
        # Apply filters button
        if st.button("Apply Filters"):
            with st.spinner("Applying filters..."):
                if SCANPY_AVAILABLE:
                    # Store original counts
                    original_cells = adata.n_obs
                    original_genes = adata.n_vars
                    
                    # Apply cell filters
                    keep_cells = True  # Default to keep all cells
                    
                    if 'n_genes_by_counts' in adata.obs.columns:
                        keep_cells = keep_cells & (adata.obs['n_genes_by_counts'] >= min_genes)
                    
                    if 'total_counts' in adata.obs.columns:
                        keep_cells = keep_cells & (adata.obs['total_counts'] >= min_counts)
                    
                    if 'pct_counts_mt' in adata.obs.columns:
                        keep_cells = keep_cells & (adata.obs['pct_counts_mt'] <= max_mt_pct)
                    
                    if 'pct_counts_hb' in adata.obs.columns:
                        keep_cells = keep_cells & (adata.obs['pct_counts_hb'] <= max_hb_pct)
                    
                    # Apply the cell filter
                    adata = adata[keep_cells, :]
                    
                    # Apply gene filters
                    keep_genes = True  # Default to keep all genes
                    
                    if filter_mt and 'mt' in adata.var.columns:
                        keep_genes = keep_genes & ~adata.var['mt']
                    
                    if filter_hb and 'hb' in adata.var.columns:
                        keep_genes = keep_genes & ~adata.var['hb']
                    
                    # Min cells per gene
                    if min_cells > 0:
                        # Calculate number of cells expressing each gene
                        n_cells_per_gene = np.sum(adata.X > 0, axis=0)
                        gene_expressed_in_min_cells = n_cells_per_gene >= (min_cells/100.0 * adata.n_obs)
                        adata.var['n_cells'] = n_cells_per_gene
                        keep_genes = keep_genes & (adata.var['n_cells'] >= min_cells/100.0 * adata.n_obs)
                    
                    # Custom genes to filter
                    if filter_custom and custom_genes:
                        gene_list = [g.strip() for g in custom_genes.split(',')]
                        keep_genes = keep_genes & ~adata.var_names.isin(gene_list)
                    
                    # Apply the gene filter
                    adata = adata[:, keep_genes]
                    
                    # Calculate how many were removed
                    cells_removed = original_cells - adata.n_obs
                    genes_removed = original_genes - adata.n_vars
                    
                    # Update the session state
                    st.session_state.data = adata
                    
                    # Save to S3
                    success, s3_uri = upload_to_s3(
                        adata,
                        f"{st.session_state.file_info['unique_name']}_filtered.h5ad",
                        QC_PREFIX
                    )
                    
                    if success:
                        st.session_state.file_info['filtered_s3_uri'] = s3_uri
                        st.success(f"Filters applied! Removed {cells_removed} cells and {genes_removed} genes.")
                        st.write(f"New dataset shape: {adata.shape}")
                    else:
                        st.error("Failed to save filtered data to S3.")
                else:
                    # Mock filtering for demonstration
                    time.sleep(2)
                    st.success("Filters applied! (Mock data)")
                    st.info("In a real application, the data would be filtered according to your criteria.")
    
    # Navigation buttons
    col1, col2 = st.columns(2)
    
    with col1:
        if st.button("⬅️ Back to Visualization"):
            st.session_state.step = 2
            st.rerun()
    
    with col2:
        if st.button("Proceed to Analysis ➡️"):
            st.session_state.step = 4
            st.rerun()

if __name__ == "__main__":
    # For testing this page in isolation
    render_qc_filtering_page()