import boto3
import os
import json
import uuid
from botocore.exceptions import NoCredentialsError
import tempfile
import io
import pickle
import streamlit as st

# S3 Prefixes for different stages
INPUT_PREFIX = "input/"
PREPARED_PREFIX = "prepared/"
VISUALIZED_PREFIX = "visualized/"
QC_PREFIX = "quality_control/"
RESULT_PREFIX = "results/"
JOB_PREFIX = "jobs/"

def init_aws_clients():
    """Initialize and return AWS client connections"""
    try:
        # Get AWS credentials from environment variables
        aws_access_key = os.getenv("AWS_ACCESS_KEY_ID")
        aws_secret_key = os.getenv("AWS_SECRET_ACCESS_KEY")
        region = os.getenv("AWS_DEFAULT_REGION", "us-east-1")
        s3_bucket = os.getenv("S3_BUCKET_NAME")
        
        # Check if credentials are available
        if not aws_access_key or not aws_secret_key or not s3_bucket:
            return {
                'error': "AWS credentials or bucket name not found in environment variables",
                'status': 'failed'
            }
        
        # Initialize clients
        s3_client = boto3.client(
            's3',
            aws_access_key_id=aws_access_key,
            aws_secret_access_key=aws_secret_key,
            region_name=region
        )
        
        sqs_client = boto3.client(
            'sqs',
            aws_access_key_id=aws_access_key,
            aws_secret_access_key=aws_secret_key,
            region_name=region
        )
        
        dynamodb_client = boto3.resource(
            'dynamodb',
            aws_access_key_id=aws_access_key,
            aws_secret_access_key=aws_secret_key,
            region_name=region
        )
        
        return {
            's3': s3_client,
            'sqs': sqs_client,
            'dynamodb': dynamodb_client,
            'bucket_name': s3_bucket
        }
        
    except Exception as e:
        # Return error details instead of raising to allow for graceful UI handling
        return {
            'error': str(e),
            'status': 'failed'
        }

def upload_to_s3(data, filename, prefix, clients=None):
    """
    Upload data to S3 with specified prefix
    
    Parameters:
    -----------
    data : object
        The data to upload (file object, AnnData, etc.)
    filename : str
        The name to give the file in S3
    prefix : str
        The S3 prefix to use
    clients : dict, optional
        Dict containing AWS clients, will use session state if not provided
        
    Returns:
    --------
    bool
        True if upload succeeded, False otherwise
    str
        S3 URI of the uploaded file if successful
    """
    try:
        import anndata as ad
    except ImportError:
        print("Warning: anndata not installed. Cannot upload AnnData objects.")
    
    # Use provided clients or get from session state
    if clients is None:
        if 'aws_clients' not in st.session_state:
            st.error("AWS clients not initialized")
            return False, None
        clients = st.session_state.aws_clients
    
    s3_client = clients['s3']
    bucket_name = clients['bucket_name']
    
    try:
        s3_key = prefix + filename
        
        # For AnnData objects, save to h5ad format
        if 'anndata' in str(type(data)).lower():
            with tempfile.NamedTemporaryFile(suffix='.h5ad') as tmp:
                data.write_h5ad(tmp.name)
                tmp.flush()  # Ensure all data is written
                tmp.seek(0)
                s3_client.upload_file(tmp.name, bucket_name, s3_key)
        
        # For file objects from file_uploader
        elif hasattr(data, 'read'):
            s3_client.upload_fileobj(data, bucket_name, s3_key)
        
        # For other data types (figures, metadata, etc.)
        else:
            buffer = io.BytesIO()
            pickle.dump(data, buffer)
            buffer.seek(0)
            s3_client.upload_fileobj(buffer, bucket_name, s3_key)
        
        # Return success and the S3 URI
        s3_uri = f"s3://{bucket_name}/{s3_key}"
        return True, s3_uri
    
    except NoCredentialsError:
        st.error("AWS credentials not found.")
        return False, None
    except Exception as e:
        st.error(f"Error uploading to S3: {e}")
        return False, None

def download_from_s3(s3_uri, clients=None):
    """
    Download data from S3
    
    Parameters:
    -----------
    s3_uri : str
        S3 URI of the file to download (s3://bucket/key)
    clients : dict, optional
        Dict containing AWS clients, will use session state if not provided
        
    Returns:
    --------
    object
        The downloaded data
    """
    try:
        import scanpy as sc
    except ImportError:
        print("Warning: scanpy not installed. Cannot download AnnData objects.")
    
    # Use provided clients or get from session state
    if clients is None:
        if 'aws_clients' not in st.session_state:
            st.error("AWS clients not initialized")
            return None
        clients = st.session_state.aws_clients
    
    s3_client = clients['s3']
    
    try:
        # Parse the S3 URI
        bucket_name = s3_uri.split('/')[2]
        s3_key = '/'.join(s3_uri.split('/')[3:])
        
        # For h5ad files, use scanpy to read
        if s3_key.endswith('.h5ad'):
            with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as tmp:
                s3_client.download_file(bucket_name, s3_key, tmp.name)
                tmp_path = tmp.name
            
            # Read with scanpy and clean up
            adata = sc.read_h5ad(tmp_path)
            os.unlink(tmp_path)
            return adata
            
        else:
            # For other files, use pickle
            buffer = io.BytesIO()
            s3_client.download_fileobj(bucket_name, s3_key, buffer)
            buffer.seek(0)
            return pickle.load(buffer)
    
    except Exception as e:
        st.error(f"Error downloading from S3: {e}")
        return None

def check_file_in_s3(s3_uri, clients=None):
    """
    Check if a file exists in S3
    
    Parameters:
    -----------
    s3_uri : str
        S3 URI of the file to check (s3://bucket/key)
    clients : dict, optional
        Dict containing AWS clients, will use session state if not provided
        
    Returns:
    --------
    bool
        True if file exists, False otherwise
    """
    # Use provided clients or get from session state
    if clients is None:
        if 'aws_clients' not in st.session_state:
            st.error("AWS clients not initialized")
            return False
        clients = st.session_state.aws_clients
    
    s3_client = clients['s3']
    
    try:
        # Parse the S3 URI
        bucket_name = s3_uri.split('/')[2]
        s3_key = '/'.join(s3_uri.split('/')[3:])
        
        # Check if file exists
        s3_client.head_object(Bucket=bucket_name, Key=s3_key)
        return True
    
    except Exception as e:
        # File doesn't exist or error occurred
        return False

def get_job_config():
    """Return the SQS queue URLs and job configurations"""
    # Get configuration from environment variables
    return {
        'queues': {
            'light_processing': os.getenv('SQS_LIGHT_QUEUE_URL'),
            'heavy_processing': os.getenv('SQS_HEAVY_QUEUE_URL')
        },
        'job_definitions': {
            'normalization': os.getenv('BATCH_NORMALIZATION_JOB_DEF'),
            'dimensionality_reduction': os.getenv('BATCH_DIM_REDUCTION_JOB_DEF'),
            'clustering': os.getenv('BATCH_CLUSTERING_JOB_DEF'),
            'differential_expression': os.getenv('BATCH_DIFF_EXP_JOB_DEF')
        },
        'dynamodb_table': os.getenv('DYNAMODB_JOB_TABLE')
    }

def get_unique_id():
    """Generate a unique ID for a job or session"""
    return str(uuid.uuid4())