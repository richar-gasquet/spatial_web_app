import json
import os
import boto3
from datetime import datetime, timedelta
from decimal import Decimal

# Initialize AWS clients
dynamodb = boto3.resource('dynamodb')
batch = boto3.client('batch')

# Get environment variables
DYNAMODB_TABLE = os.environ.get('DYNAMODB_TABLE')

# DynamoDB doesn't support float/int types directly, we need to convert them to decimal
class DecimalEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Decimal):
            return float(obj)
        return super(DecimalEncoder, self).default(obj)

def lambda_handler(event, context):
    """
    Lambda handler for periodically checking job status
    
    Parameters:
    -----------
    event : dict
        The CloudWatch scheduled event trigger
    context : LambdaContext
        The Lambda context
        
    Returns:
    --------
    dict
        Response indicating success or failure
    """
    # Get the DynamoDB table
    table = dynamodb.Table(DYNAMODB_TABLE)
    
    # Get all jobs that are in progress (submitted to batch but not complete)
    response = table.scan(
        FilterExpression="(#status = :submitted OR #status = :running) AND attribute_exists(batch_job_id)",
        ExpressionAttributeNames={'#status': 'status'},
        ExpressionAttributeValues={
            ':submitted': 'SUBMITTED_TO_BATCH', 
            ':running': 'RUNNING'
        }
    )
    
    # Process each in-progress job
    jobs_processed = 0
    jobs_updated = 0
    
    for job in response.get('Items', []):
        jobs_processed += 1
        job_id = job.get('job_id')
        batch_job_id = job.get('batch_job_id')
        
        if not batch_job_id:
            print(f"Job {job_id} has no batch_job_id, skipping")
            continue
        
        try:
            # Get the batch job status
            batch_response = batch.describe_jobs(jobs=[batch_job_id])
            
            # If job not found in Batch, mark as failed
            if not batch_response.get('jobs'):
                print(f"Job {job_id} (Batch ID: {batch_job_id}) not found in Batch")
                table.update_item(
                    Key={'job_id': job_id},
                    UpdateExpression="SET #status = :status, updated_at = :updated_at, error = :error",
                    ExpressionAttributeNames={'#status': 'status'},
                    ExpressionAttributeValues={
                        ':status': 'FAILED',
                        ':updated_at': datetime.now().isoformat(),
                        ':error': 'Job not found in AWS Batch'
                    }
                )
                jobs_updated += 1
                continue
            
            # Get the first job (should only be one)
            batch_job = batch_response['jobs'][0]
            batch_status = batch_job['status']
            
            # Map Batch status to our status
            status_mapping = {
                'SUBMITTED': 'SUBMITTED_TO_BATCH',
                'PENDING': 'SUBMITTED_TO_BATCH',
                'RUNNABLE': 'SUBMITTED_TO_BATCH',
                'STARTING': 'RUNNING',
                'RUNNING': 'RUNNING',
                'SUCCEEDED': 'SUCCEEDED',
                'FAILED': 'FAILED'
            }
            
            # Get our status from mapping
            new_status = status_mapping.get(batch_status, 'UNKNOWN')
            
            # Update the job status if it changed
            if new_status != job.get('status'):
                update_expr = "SET #status = :status, updated_at = :updated_at"
                expr_values = {
                    ':status': new_status,
                    ':updated_at': datetime.now().isoformat()
                }
                
                # If the job failed, include the reason
                if new_status == 'FAILED':
                    update_expr += ", error = :error"
                    container = batch_job.get('container', {})
                    reason = container.get('reason', 'Unknown error')
                    expr_values[':error'] = reason
                
                # If the job succeeded, include the result S3 URI if available
                if new_status == 'SUCCEEDED':
                    # Check for environment variables that might contain result information
                    for env in batch_job.get('container', {}).get('environment', []):
                        if env['name'] == 'OUTPUT_S3_URI':
                            update_expr += ", result_s3_uri = :result_uri"
                            expr_values[':result_uri'] = env['value']
                    
                    # Add completed time
                    update_expr += ", completed_at = :completed_at"
                    expr_values[':completed_at'] = datetime.now().isoformat()
                
                # Update the item in DynamoDB
                table.update_item(
                    Key={'job_id': job_id},
                    UpdateExpression=update_expr,
                    ExpressionAttributeNames={'#status': 'status'},
                    ExpressionAttributeValues=expr_values
                )
                
                jobs_updated += 1
                print(f"Updated job {job_id} status from {job.get('status')} to {new_status}")
            
        except Exception as e:
            print(f"Error processing job {job_id}: {e}")
    
    # Clean up old jobs (optional)
    try:
        # Find jobs older than 7 days
        cutoff_date = (datetime.now() - timedelta(days=7)).isoformat()
        
        old_jobs_response = table.scan(
            FilterExpression="(#status = :succeeded OR #status = :failed) AND updated_at < :cutoff",
            ExpressionAttributeNames={'#status': 'status'},
            ExpressionAttributeValues={
                ':succeeded': 'SUCCEEDED',
                ':failed': 'FAILED',
                ':cutoff': cutoff_date
            }
        )
        
        # Just log old jobs for now (could delete them if needed)
        old_jobs = len(old_jobs_response.get('Items', []))
        if old_jobs > 0:
            print(f"Found {old_jobs} completed jobs older than 7 days")
    
    except Exception as e:
        print(f"Error during cleanup: {e}")
    
    return {
        'statusCode': 200,
        'body': json.dumps({
            'jobs_processed': jobs_processed,
            'jobs_updated': jobs_updated
        }, cls=DecimalEncoder)
    }