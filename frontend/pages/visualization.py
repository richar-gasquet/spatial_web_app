import streamlit as st
import pandas as pd
import numpy as np
import io
import time

from utils.aws_utils import download_from_s3, upload_to_s3, VISUALIZED_PREFIX
from utils.ui_components import render_data_info

# Try to import visualization libraries
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scanpy as sc
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False
    print("Warning: visualization libraries not available.")

def render_visualization_page():
    """Render the data visualization page"""
    st.header("📈 Data Visualization")
    
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
    st.markdown("### Visualization Options")
    
    # Visualization type selector
    viz_type = st.selectbox(
        "Visualization Type",
        ["QC Metrics", "Gene Expression", "Spatial Distribution"]
    )
    
    # Different visualizations based on selection
    if viz_type == "QC Metrics":
        render_qc_visualizations(adata)
    elif viz_type == "Gene Expression":
        render_gene_visualizations(adata)
    elif viz_type == "Spatial Distribution":
        render_spatial_visualizations(adata)
    
    # Save visualization button
    if st.button("Save Visualization Results"):
        with st.spinner("Saving visualization results..."):
            # Save the dataset with visualization results to S3
            if PLOTTING_AVAILABLE:
                success, s3_uri = upload_to_s3(
                    adata,
                    f"{st.session_state.file_info['unique_name']}_visualized.h5ad",
                    VISUALIZED_PREFIX
                )
                
                if success:
                    st.success("Visualization results saved to S3!")
                    st.session_state.file_info['visualized_s3_uri'] = s3_uri
                else:
                    st.error("Failed to save visualization results to S3.")
            else:
                # Mock saving for demonstration
                time.sleep(1)
                st.success("Visualization results saved (mock).")
                st.session_state.file_info['visualized_s3_uri'] = "mock_uri"
    
    # Navigation buttons
    col1, col2 = st.columns(2)
    
    with col1:
        if st.button("⬅️ Back to Preparation"):
            st.session_state.step = 1
            st.rerun()
    
    with col2:
        if st.button("Proceed to Quality Control ➡️"):
            st.session_state.step = 3
            st.rerun()

def render_qc_visualizations(adata):
    """Render QC metric visualizations"""
    if not hasattr(adata, 'obs') or 'n_genes_by_counts' not in adata.obs.columns:
        st.warning("QC metrics not found. Please go back and calculate QC metrics.")
        return
    
    # Get available QC metrics
    qc_cols = [col for col in adata.obs.columns if any(x in col for x in ['counts', 'pct', 'n_genes'])]
    
    if not qc_cols:
        st.warning("No QC metrics found in the dataset.")
        return
    
    # Let user select metrics to visualize
    selected_metrics = st.multiselect(
        "Select QC Metrics to Visualize",
        qc_cols,
        default=qc_cols[:3] if len(qc_cols) >= 3 else qc_cols
    )
    
    if not selected_metrics:
        st.info("Please select at least one metric to visualize.")
        return
    
    # Display visualizations
    if PLOTTING_AVAILABLE:
        # Create violin plots
        try:
            # Using scanpy's plotting functions
            for metric in selected_metrics:
                st.subheader(f"{metric} Distribution")
                fig = plt.figure(figsize=(10, 6))
                if hasattr(adata.obs, 'library_id'):
                    sc.pl.violin(adata, metric, groupby='library_id', ax=fig.gca(), show=False)
                else:
                    sc.pl.violin(adata, metric, ax=fig.gca(), show=False)
                st.pyplot(fig)
                plt.close(fig)
            
            # Scatter plots for relationships between metrics
            if len(selected_metrics) >= 2:
                st.subheader("Relationships Between Metrics")
                fig = plt.figure(figsize=(10, 6))
                sns.scatterplot(
                    data=adata.obs, 
                    x=selected_metrics[0], 
                    y=selected_metrics[1],
                    alpha=0.5
                )
                plt.title(f"{selected_metrics[0]} vs {selected_metrics[1]}")
                st.pyplot(fig)
                plt.close(fig)
        except Exception as e:
            st.error(f"Error creating QC plots: {e}")
    else:
        # Show mock visualizations
        st.info("QC visualizations would be shown here. Install visualization libraries to see actual plots.")
        
        # Show data table instead
        st.write("Sample of QC metrics:")
        st.dataframe(adata.obs[selected_metrics].head(10))

def render_gene_visualizations(adata):
    """Render gene expression visualizations"""
    # Get available genes
    if hasattr(adata, 'var_names'):
        all_genes = adata.var_names.tolist()
        # Show only first 1000 genes if there are too many
        if len(all_genes) > 1000:
            st.info(f"Showing first 1000 of {len(all_genes)} genes.")
            gene_list = all_genes[:1000]
        else:
            gene_list = all_genes
    else:
        st.warning("No gene names found in the dataset.")
        return
    
    # Let user search for genes
    gene_search = st.text_input("Search for genes", "")
    
    if gene_search:
        matching_genes = [gene for gene in all_genes if gene_search.lower() in gene.lower()]
        if not matching_genes:
            st.warning(f"No genes found matching '{gene_search}'")
        else:
            st.info(f"Found {len(matching_genes)} matching genes.")
            gene_list = matching_genes[:1000]  # Limit to 1000 matches
    
    # Let user select genes to visualize
    selected_genes = st.multiselect(
        "Select Genes to Visualize",
        gene_list,
        default=gene_list[:3] if len(gene_list) >= 3 else gene_list
    )
    
    if not selected_genes:
        st.info("Please select at least one gene to visualize.")
        return
    
    # Display visualizations
    if PLOTTING_AVAILABLE:
        try:
            # Display expression levels
            for gene in selected_genes:
                st.subheader(f"{gene} Expression")
                
                # Check if spatial coordinates are available
                if 'spatial' in adata.obsm:
                    # Spatial plot
                    fig = plt.figure(figsize=(10, 10))
                    
                    if adata.file_info['file_type'] == 'visium' and 'spatial' in adata.uns:
                        # Visium spatial plot
                        try:
                            sc.pl.spatial(adata, color=gene, show=False, ax=fig.gca())
                        except Exception as e:
                            st.error(f"Error creating spatial plot: {e}")
                    else:
                        # Generic spatial scatter plot
                        coords = adata.obsm['spatial']
                        plt.scatter(
                            coords[:, 0], 
                            coords[:, 1], 
                            c=adata[:, gene].X.flatten(),
                            cmap='viridis',
                            s=10
                        )
                        plt.colorbar(label='Expression')
                        plt.title(f"{gene} Expression")
                    
                    st.pyplot(fig)
                    plt.close(fig)
                else:
                    # Violin plot as fallback
                    fig = plt.figure(figsize=(10, 6))
                    sc.pl.violin(adata, [gene], ax=fig.gca(), show=False)
                    st.pyplot(fig)
                    plt.close(fig)
        except Exception as e:
            st.error(f"Error creating gene expression plots: {e}")
    else:
        # Show mock visualizations
        st.info("Gene expression visualizations would be shown here. Install visualization libraries to see actual plots.")

def render_spatial_visualizations(adata):
    """Render spatial distribution visualizations"""
    # Check if spatial coordinates are available
    if not hasattr(adata, 'obsm') or 'spatial' not in adata.obsm:
        st.warning("No spatial coordinates found in the dataset.")
        return
    
    # Let user select visualization options
    st.subheader("Spatial Distribution Options")
    
    # Get categorical columns for coloring
    cat_cols = [
        col for col in adata.obs.columns 
        if adata.obs[col].dtype.name == 'category' or 
        len(adata.obs[col].unique()) < 20
    ]
    
    color_by = st.selectbox(
        "Color by",
        ["None"] + cat_cols,
        index=0
    )
    
    # Display visualizations
    if PLOTTING_AVAILABLE:
        try:
            st.subheader("Spatial Distribution")
            
            # Check if it's a Visium dataset with image
            if adata.file_info['file_type'] == 'visium' and 'spatial' in adata.uns:
                # Visium spatial plot with tissue image
                fig = plt.figure(figsize=(12, 12))
                
                if color_by != "None":
                    sc.pl.spatial(adata, color=color_by, show=False, ax=fig.gca())
                else:
                    sc.pl.spatial(adata, show=False, ax=fig.gca())
                
                st.pyplot(fig)
                plt.close(fig)
            else:
                # Generic spatial scatter plot
                fig = plt.figure(figsize=(10, 10))
                
                coords = adata.obsm['spatial']
                
                if color_by != "None":
                    if adata.obs[color_by].dtype.name == 'category':
                        # Categorical coloring
                        categories = adata.obs[color_by].cat.categories
                        cmap = plt.cm.get_cmap('tab20', len(categories))
                        
                        for i, cat in enumerate(categories):
                            mask = adata.obs[color_by] == cat
                            plt.scatter(
                                coords[mask, 0], 
                                coords[mask, 1],
                                c=[cmap(i)],
                                label=cat,
                                alpha=0.7,
                                s=10
                            )
                        plt.legend(title=color_by)
                    else:
                        # Continuous coloring
                        plt.scatter(
                            coords[:, 0], 
                            coords[:, 1],
                            c=adata.obs[color_by],
                            cmap='viridis',
                            alpha=0.7,
                            s=10
                        )
                        plt.colorbar(label=color_by)
                else:
                    # No coloring
                    plt.scatter(
                        coords[:, 0], 
                        coords[:, 1],
                        alpha=0.7,
                        s=10
                    )
                
                plt.title("Spatial Distribution")
                plt.axis('equal')
                st.pyplot(fig)
                plt.close(fig)
        except Exception as e:
            st.error(f"Error creating spatial plots: {e}")
    else:
        # Show mock visualizations
        st.info("Spatial visualizations would be shown here. Install visualization libraries to see actual plots.")

if __name__ == "__main__":
    # For testing this page in isolation
    render_visualization_page()