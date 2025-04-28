import json
import os
import boto3
import uuid
from datetime import datetime

# Initialize AWS clients
dynamodb = boto3.resource('dynamodb')
s3 = boto3.client('s3')
batch = boto3.client('batch')

# Get environment variables
DYNAMODB_TABLE = os.environ.get('DYNAMODB_TABLE')
BATCH_JOB_QUEUE = os.environ.get('BATCH_JOB_QUEUE')
BATCH_JOB_DEFINITION_PREFIX = os.environ.get('BATCH_JOB_DEFINITION_PREFIX', 'spatial-transcriptomics-')

# Define job definition mapping
JOB_DEFINITION_MAPPING = {
    'normalization': 'normalization',
    'normalization_hvg_selection': 'normalization',
    'dimensionality_reduction': 'dimensionality_reduction',
    'pca_umap': 'dimensionality_reduction',
    'clustering': 'clustering',
    'differential_expression': 'differential_expression'
}

def lambda_handler(event, context):
    """
    Lambda handler for processing SQS messages and submitting Batch jobs
    
    Parameters:
    -----------
    event : dict
        The event dict containing the SQS message
    context : LambdaContext
        The Lambda context
        
    Returns:
    --------
    dict
        Response indicating success or failure
    """
    print(f"Received event: {json.dumps(event)}")
    
    # Process SQS messages
    for record in event.get('Records', []):
        try:
            # Parse the message body
            message = json.loads(record['body'])
            
            # Extract job details
            job_id = message.get('job_id')
            job_type = message.get('job_type')
            input_data = message.get('input_data', {})
            job_params = message.get('job_params', {})
            
            # Validate required fields
            if not job_id or not job_type or not input_data:
                print(f"Invalid message: missing required fields: {message}")
                continue
            
            # Get job definition name
            job_definition_key = JOB_DEFINITION_MAPPING.get(job_type, job_type)
            job_definition = f"{BATCH_JOB_DEFINITION_PREFIX}{job_definition_key}"
            
            # Update job status in DynamoDB
            table = dynamodb.Table(DYNAMODB_TABLE)
            
            # Mark job as submitted to batch
            update_result = table.update_item(
                Key={'job_id': job_id},
                UpdateExpression="SET #status = :status, updated_at = :updated_at",
                ExpressionAttributeNames={'#status': 'status'},
                ExpressionAttributeValues={
                    ':status': 'SUBMITTED_TO_BATCH',
                    ':updated_at': datetime.now().isoformat()
                },
                ReturnValues="UPDATED_NEW"
            )
            
            print(f"Updated job status in DynamoDB: {update_result}")
            
            # Prepare container override environment variables
            environment = [
                {'name': 'JOB_ID', 'value': job_id},
                {'name': 'JOB_TYPE', 'value': job_type},
                {'name': 'DYNAMODB_TABLE', 'value': DYNAMODB_TABLE},
                {'name': 'INPUT_S3_URI', 'value': input_data.get('s3_uri', '')},
                {'name': 'OUTPUT_FILENAME', 'value': input_data.get('filename', 'output')}
            ]
            
            # Add job parameters as environment variables
            for key, value in job_params.items():
                if isinstance(value, (str, int, float, bool)):
                    environment.append({
                        'name': f"PARAM_{key.upper()}",
                        'value': str(value)
                    })
                else:
                    # For complex types (lists, dicts), convert to JSON
                    environment.append({
                        'name': f"PARAM_{key.upper()}",
                        'value': json.dumps(value)
                    })
            
            # Submit the job to AWS Batch
            response = batch.submit_job(
                jobName=f"{job_type}-{job_id}",
                jobQueue=BATCH_JOB_QUEUE,
                jobDefinition=job_definition,
                containerOverrides={
                    'environment': environment
                }
            )
            
            # Get the AWS Batch job ID
            batch_job_id = response['jobId']
            
            # Update job in DynamoDB with Batch job ID
            table.update_item(
                Key={'job_id': job_id},
                UpdateExpression="SET batch_job_id = :batch_job_id",
                ExpressionAttributeValues={
                    ':batch_job_id': batch_job_id
                }
            )
            
            print(f"Submitted job to AWS Batch: {batch_job_id}")
            
        except Exception as e:
            print(f"Error processing message: {e}")
            continue
    
    return {
        'statusCode': 200,
        'body': json.dumps('Job processing complete')
    }