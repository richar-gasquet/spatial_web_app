import streamlit as st
import os
import sys
import json
from datetime import datetime

# Add the project root to the path so we can import utils
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import utility modules
from utils.aws_utils import init_aws_clients
from utils.ui_components import render_sidebar, render_header

# Load environment variables
from dotenv import load_dotenv
load_dotenv()

# Set page configuration
st.set_page_config(
    page_title="Spatial Transcriptomics Analysis",
    page_icon="📊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Initialize session state
if 'initialized' not in st.session_state:
    st.session_state.initialized = True
    st.session_state.step = 0
    st.session_state.file_info = {
        'original_name': None,
        'unique_name': None,
        'file_type': None,
        'file_path': None
    }
    st.session_state.jobs = {}  # Store job information
    st.session_state.data = None  # For storing the current dataset
    
    # Initialize AWS clients
    st.session_state.aws_clients = init_aws_clients()

# Application title and header
render_header()

# Sidebar navigation
render_sidebar()

# Main content based on current step
steps = [
    "Upload Data",
    "Data Preparation",
    "Data Visualization", 
    "Quality Control & Filtering",
    "Analysis & Results"
]

# Current step content
if st.session_state.step == 0:
    # Upload step
    from pages.upload import render_upload_page
    render_upload_page()
    
elif st.session_state.step == 1:
    # Preparation step
    from pages.preparation import render_preparation_page
    render_preparation_page()
    
elif st.session_state.step == 2:
    # Visualization step
    from pages.visualization import render_visualization_page
    render_visualization_page()
    
elif st.session_state.step == 3:
    # QC & Filtering step
    from pages.qc_filtering import render_qc_filtering_page
    render_qc_filtering_page()
    
elif st.session_state.step == 4:
    # Analysis step
    from pages.analysis import render_analysis_page
    render_analysis_page()

# App footer
st.markdown("---")
st.markdown("*Spatial Transcriptomics Analysis Pipeline* | Developed with Streamlit & Scanpy")