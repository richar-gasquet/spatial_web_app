import os
import json
import boto3
import logging
from datetime import datetime
import tempfile
import scanpy as sc
import anndata as ad

# Configure logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('batch_utils')

# Initialize AWS clients
s3 = boto3.client('s3')
dynamodb = boto3.resource('dynamodb')

def update_job_status(status, additional_data=None):
    """
    Update job status in DynamoDB
    
    Parameters:
    -----------
    status : str
        New status (RUNNING, SUCCEEDED, FAILED)
    additional_data : dict, optional
        Additional data to store with the job
    """
    job_id = os.environ.get('JOB_ID')
    table_name = os.environ.get('DYNAMODB_TABLE')
    
    if not job_id or not table_name:
        logger.error("Missing JOB_ID or DYNAMODB_TABLE environment variables")
        return False
    
    try:
        table = dynamodb.Table(table_name)
        
        update_expr = "SET #status = :status, updated_at = :updated_at"
        expr_values = {
            ':status': status,
            ':updated_at': datetime.now().isoformat()
        }
        
        # Add additional data if provided
        if additional_data and isinstance(additional_data, dict):
            for key, value in additional_data.items():
                update_expr += f", {key} = :{key}"
                expr_values[f":{key}"] = value
        
        response = table.update_item(
            Key={'job_id': job_id},
            UpdateExpression=update_expr,
            ExpressionAttributeNames={'#status': 'status'},
            ExpressionAttributeValues=expr_values,
            ReturnValues="UPDATED_NEW"
        )
        
        logger.info(f"Updated job status to {status}")
        return True
    
    except Exception as e:
        logger.error(f"Error updating job status: {e}")
        return False

def get_input_data():
    """
    Download input data from S3
    
    Returns:
    --------
    anndata.AnnData
        The loaded data object
    """
    input_s3_uri = os.environ.get('INPUT_S3_URI')
    
    if not input_s3_uri:
        logger.error("Missing INPUT_S3_URI environment variable")
        raise ValueError("No input S3 URI provided")
    
    # Parse S3 URI
    if not input_s3_uri.startswith('s3://'):
        raise ValueError(f"Invalid S3 URI: {input_s3_uri}")
    
    parts = input_s3_uri.replace('s3://', '').split('/')
    bucket = parts[0]
    key = '/'.join(parts[1:])
    
    logger.info(f"Downloading from S3: {bucket}/{key}")
    
    # Download to temporary file
    with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as tmp:
        s3.download_file(bucket, key, tmp.name)
        tmp_path = tmp.name
    
    # Load with scanpy
    try:
        adata = sc.read_h5ad(tmp_path)
        os.unlink(tmp_path)  # Delete the temp file
        logger.info(f"Successfully loaded data: {adata.shape}")
        return adata
    except Exception as e:
        logger.error(f"Error loading data: {e}")
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise

def save_output_data(adata, additional_outputs=None):
    """
    Save output data to S3
    
    Parameters:
    -----------
    adata : anndata.AnnData
        The data to save
    additional_outputs : dict, optional
        Additional outputs to save (e.g., results, plots)
        
    Returns:
    --------
    dict
        Dictionary with S3 URIs of saved objects
    """
    output_filename = os.environ.get('OUTPUT_FILENAME', 'output')
    job_id = os.environ.get('JOB_ID', 'unknown')
    
    try:
        # Determine bucket from input URI
        input_s3_uri = os.environ.get('INPUT_S3_URI')
        bucket = input_s3_uri.replace('s3://', '').split('/')[0]
        
        # Save to temporary file then upload
        with tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False) as tmp:
            adata.write_h5ad(tmp.name)
            tmp_path = tmp.name
        
        # Construct output key
        output_key = f"results/{output_filename}_{job_id}.h5ad"
        
        # Upload to S3
        s3.upload_file(tmp_path, bucket, output_key)
        os.unlink(tmp_path)
        
        output_s3_uri = f"s3://{bucket}/{output_key}"
        logger.info(f"Saved output to S3: {output_s3_uri}")
        
        result = {'output_s3_uri': output_s3_uri}
        
        # Save additional outputs if provided
        if additional_outputs and isinstance(additional_outputs, dict):
            for name, data in additional_outputs.items():
                if isinstance(data, (dict, list)):
                    # Save JSON data
                    tmp_output_key = f"results/{output_filename}_{job_id}_{name}.json"
                    s3.put_object(
                        Body=json.dumps(data),
                        Bucket=bucket,
                        Key=tmp_output_key
                    )
                    result[f"{name}_s3_uri"] = f"s3://{bucket}/{tmp_output_key}"
                else:
                    # For other types, use pickle (already implemented in frontend utils)
                    logger.warning(f"Skipping additional output {name} - unsupported type")
        
        # Set environment variable for status updater to read
        os.environ['OUTPUT_S3_URI'] = output_s3_uri
        
        return result
        
    except Exception as e:
        logger.error(f"Error saving output: {e}")
        raise

def get_job_params():
    """
    Get job parameters from environment variables
    
    Returns:
    --------
    dict
        Dictionary of job parameters
    """
    params = {}
    
    # Get all environment variables starting with PARAM_
    for key, value in os.environ.items():
        if key.startswith('PARAM_'):
            param_name = key[6:].lower()  # Remove PARAM_ and convert to lowercase
            
            # Try to parse JSON for complex types
            try:
                params[param_name] = json.loads(value)
            except json.JSONDecodeError:
                # Handle primitive types
                if value.lower() == 'true':
                    params[param_name] = True
                elif value.lower() == 'false':
                    params[param_name] = False
                elif value.isdigit():
                    params[param_name] = int(value)
                elif value.replace('.', '').isdigit() and value.count('.') == 1:
                    params[param_name] = float(value)
                else:
                    params[param_name] = value
    
    return params