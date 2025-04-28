import streamlit as st
import os
import tempfile
import time
from utils.aws_utils import upload_to_s3, INPUT_PREFIX, get_unique_id

def render_upload_page():
    """Render the data upload page"""
    st.header("📤 Upload Data")
    
    st.markdown("""
    Upload your spatial transcriptomics data file. Supported formats:
    - 10x Genomics Visium data (folder containing h5 matrices)
    - H5AD files (pre-processed data)
    - CSV/TSV matrices
    """)
    
    # File type selection
    col1, col2 = st.columns(2)
    with col1:
        file_type = st.selectbox(
            "Data Type",
            ["visium", "h5ad", "csv/tsv"],
            help="Select the type of data you are uploading"
        )
    
    # Optional sample name
    with col2:
        sample_name = st.text_input(
            "Sample Name (optional)",
            help="Give your sample a descriptive name"
        )
    
    # File uploader
    uploaded_file = st.file_uploader("Choose a file to upload", type=None)
    
    if uploaded_file is not None:
        # Display file info
        st.write(f"File: `{uploaded_file.name}` ({uploaded_file.size / 1e6:.2f} MB)")
        
        # Generate a unique filename
        if sample_name:
            unique_filename = f"{sample_name}_{get_unique_id()}"
        else:
            unique_filename = f"{uploaded_file.name}_{get_unique_id()}"
        
        # Update session state with file info
        st.session_state.file_info = {
            'original_name': uploaded_file.name,
            'unique_name': unique_filename,
            'file_type': file_type,
            'file_size_mb': uploaded_file.size / 1e6
        }
        
        # Display upload status
        with st.spinner(f"Uploading `{uploaded_file.name}` to S3..."):
            # Upload to S3
            success, s3_uri = upload_to_s3(uploaded_file, unique_filename, INPUT_PREFIX)
            
            if success:
                st.success("File uploaded successfully!")
                st.session_state.file_info['s3_uri'] = s3_uri
                
                # Loading indicator
                with st.spinner("Initializing dataset..."):
                    # This is where we would load the data for validation
                    # For this example, we'll simulate a delay
                    time.sleep(1)
                    
                    # Mock successful loading
                    st.success("Dataset initialized successfully!")
                    
                    # Enable proceeding to next step
                    if st.button("Proceed to Data Preparation"):
                        st.session_state.step = 1
                        st.rerun()
            else:
                st.error("Failed to upload file to S3. Please check your AWS credentials.")
    
    # Show demo data option
    st.markdown("---")
    st.markdown("### 📊 Or use demo data")
    
    if st.button("Load Demo Dataset"):
        # Mock loading demo data
        with st.spinner("Loading demo dataset..."):
            # Simulate delay
            time.sleep(1)
            
            # Set file info for demo data
            st.session_state.file_info = {
                'original_name': "visium_mouse_brain_demo.h5ad",
                'unique_name': f"demo_mouse_brain_{get_unique_id()}",
                'file_type': "visium",
                'file_size_mb': 120.5,
                's3_uri': "s3://demo-bucket/input/visium_mouse_brain_demo.h5ad"  # Mock URI
            }
            
            st.success("Demo dataset loaded successfully!")
            
            # Enable proceeding to next step
            if st.button("Proceed to Data Preparation", key="demo_proceed"):
                st.session_state.step = 1
                st.rerun()

if __name__ == "__main__":
    # For testing this page in isolation
    render_upload_page()