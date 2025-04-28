import json
import time
import os
from datetime import datetime
import streamlit as st

from utils.aws_utils import get_job_config, get_unique_id, download_from_s3, upload_to_s3

def determine_job_type(data_size, task_type):
    """
    Determine if a job should be processed locally or in AWS Batch
    
    Parameters:
    -----------
    data_size : int
        Size of the data in MB
    task_type : str
        Type of task (normalization, clustering, etc.)
        
    Returns:
    --------
    str
        'local', 'queue_light', or 'queue_heavy'
    """
    # These thresholds would be determined based on testing and resource considerations
    
    # Define computation-intensive tasks
    heavy_tasks = [
        'dimensionality_reduction',
        'clustering',
        'differential_expression'
    ]
    
    # For very small data, always process locally
    if data_size < 5:  # Less than 5MB
        return 'local'
    
    # For medium data, use light queue unless it's a heavy task
    if data_size < 50:  # 5-50MB
        if task_type in heavy_tasks:
            return 'queue_light'
        else:
            return 'local'
    
    # For large data, use heavy queue for heavy tasks, light queue for others
    if task_type in heavy_tasks:
        return 'queue_heavy'
    else:
        return 'queue_light'

def submit_job(job_type, input_data, job_params):
    """
    Submit a job for processing
    
    Parameters:
    -----------
    job_type : str
        Type of job (e.g., 'normalization', 'clustering')
    input_data : dict
        Dict containing input data info (s3_uri, etc.)
    job_params : dict
        Parameters for the job
        
    Returns:
    --------
    dict
        Job submission details including job_id
    """
    # Get configuration
    config = get_job_config()
    
    # Generate a unique job ID
    job_id = get_unique_id()
    
    # Create the job message
    job_message = {
        'job_id': job_id,
        'job_type': job_type,
        'input_data': input_data,
        'job_params': job_params,
        'submitted_at': datetime.now().isoformat(),
        'status': 'SUBMITTED'
    }
    
    # Determine which queue to use
    processing_type = determine_job_type(
        input_data.get('data_size_mb', 100),
        job_type
    )
    
    # Submit to appropriate queue or process locally
    if processing_type == 'local':
        # For local processing, we'll add to session state and process in the Streamlit app
        if 'local_jobs' not in st.session_state:
            st.session_state.local_jobs = {}
        
        st.session_state.local_jobs[job_id] = job_message
        return {
            'job_id': job_id,
            'processing_type': 'local',
            'status': 'SUBMITTED'
        }
    else:
        # For queue processing, submit to SQS
        queue_url = config['queues']['light_processing'] if processing_type == 'queue_light' else config['queues']['heavy_processing']
        
        # Get SQS client
        sqs_client = st.session_state.aws_clients['sqs']
        
        # Send to SQS
        try:
            response = sqs_client.send_message(
                QueueUrl=queue_url,
                MessageBody=json.dumps(job_message)
            )
            
            # Also store the job in DynamoDB for tracking
            dynamodb = st.session_state.aws_clients['dynamodb']
            table = dynamodb.Table(config['dynamodb_table'])
            
            table.put_item(
                Item=job_message
            )
            
            return {
                'job_id': job_id,
                'processing_type': processing_type,
                'sqs_message_id': response.get('MessageId'),
                'status': 'SUBMITTED'
            }
        except Exception as e:
            st.error(f"Error submitting job to queue: {e}")
            return {
                'error': str(e),
                'status': 'FAILED'
            }

def check_job_status(job_id):
    """
    Check the status of a job
    
    Parameters:
    -----------
    job_id : str
        The ID of the job to check
        
    Returns:
    --------
    dict
        Job status details
    """
    # For local jobs, check session state
    if 'local_jobs' in st.session_state and job_id in st.session_state.local_jobs:
        return st.session_state.local_jobs[job_id]
    
    # For queue jobs, check DynamoDB
    config = get_job_config()
    dynamodb = st.session_state.aws_clients['dynamodb']
    table = dynamodb.Table(config['dynamodb_table'])
    
    try:
        response = table.get_item(
            Key={'job_id': job_id}
        )
        
        if 'Item' in response:
            return response['Item']
    except Exception as e:
        st.error(f"Error checking job status: {e}")
    
    return {'job_id': job_id, 'status': 'NOT_FOUND'}

def process_local_job(job_id):
    """
    Process a job locally within the Streamlit app
    
    Parameters:
    -----------
    job_id : str
        The ID of the job to process
        
    Returns:
    --------
    dict
        Job result
    """
    if 'local_jobs' not in st.session_state or job_id not in st.session_state.local_jobs:
        return {'status': 'NOT_FOUND'}
    
    job = st.session_state.local_jobs[job_id]
    
    # Update status
    job['status'] = 'RUNNING'
    st.session_state.local_jobs[job_id] = job
    
    try:
        # Depending on job_type, call appropriate processing function
        job_type = job['job_type']
        input_data = job['input_data']
        job_params = job['job_params']
        
        # Import processing functions
        from utils.data_processing import (
            process_normalization,
            process_dimensionality_reduction,
            process_clustering,
            process_differential_expression
        )
        
        # Call appropriate function
        if job_type == 'normalization':
            result = process_normalization(input_data, job_params)
        elif job_type == 'dimensionality_reduction':
            result = process_dimensionality_reduction(input_data, job_params)
        elif job_type == 'clustering':
            result = process_clustering(input_data, job_params)
        elif job_type == 'differential_expression':
            result = process_differential_expression(input_data, job_params)
        else:
            raise ValueError(f"Unknown job type: {job_type}")
        
        # Update job with result
        job['status'] = 'SUCCEEDED'
        job['completed_at'] = datetime.now().isoformat()
        job['result'] = result
        st.session_state.local_jobs[job_id] = job
        
        return job
    
    except Exception as e:
        # Update job with error
        job['status'] = 'FAILED'
        job['error'] = str(e)
        job['completed_at'] = datetime.now().isoformat()
        st.session_state.local_jobs[job_id] = job
        
        return job

def get_job_results(job_id):
    """
    Get the results of a completed job
    
    Parameters:
    -----------
    job_id : str
        The ID of the job
        
    Returns:
    --------
    dict
        Job results
    """
    # Check job status
    job_status = check_job_status(job_id)
    
    if job_status['status'] != 'SUCCEEDED':
        return {
            'status': job_status['status'],
            'error': job_status.get('error', 'Job not complete')
        }
    
    # For local jobs, get from session state
    if 'local_jobs' in st.session_state and job_id in st.session_state.local_jobs:
        return st.session_state.local_jobs[job_id]['result']
    
    # For queue jobs, get result from S3
    result_s3_uri = job_status.get('result_s3_uri')
    
    if not result_s3_uri:
        return {'status': 'ERROR', 'error': 'No result found for job'}
    
    # Download result from S3
    result = download_from_s3(result_s3_uri)
    
    return result