# AWS Deployment Instructions

This document provides step-by-step instructions for setting up the AWS infrastructure required for the backend processing of the Spatial Transcriptomics Analysis Pipeline.

## Prerequisites

1. AWS CLI installed and configured with appropriate credentials
2. Docker installed locally
3. Access to create AWS resources (S3, SQS, Lambda, DynamoDB, Batch)

## Step 1: Create S3 Bucket

First, create an S3 bucket to store your data:

```bash
aws s3api create-bucket \
    --bucket spatial-transcriptomics-data \
    --region us-east-1
```

Enable versioning on the bucket:

```bash
aws s3api put-bucket-versioning \
    --bucket spatial-transcriptomics-data \
    --versioning-configuration Status=Enabled
```

## Step 2: Create DynamoDB Table

Create a DynamoDB table for job tracking:

```bash
aws dynamodb create-table \
    --table-name spatial-transcriptomics-jobs \
    --attribute-definitions AttributeName=job_id,AttributeType=S \
    --key-schema AttributeName=job_id,KeyType=HASH \
    --billing-mode PAY_PER_REQUEST
```

## Step 3: Create SQS Queues

Create SQS queues for job processing:

```bash
# Create light processing queue
aws sqs create-queue \
    --queue-name spatial-transcriptomics-light-processing

# Create heavy processing queue
aws sqs create-queue \
    --queue-name spatial-transcriptomics-heavy-processing
```

Get the queue URLs:

```bash
aws sqs get-queue-url --queue-name spatial-transcriptomics-light-processing
aws sqs get-queue-url --queue-name spatial-transcriptomics-heavy-processing
```

Add these queue URLs to your `.env` file:

```
SQS_LIGHT_QUEUE_URL=https://sqs.region.amazonaws.com/account-id/spatial-transcriptomics-light-processing
SQS_HEAVY_QUEUE_URL=https://sqs.region.amazonaws.com/account-id/spatial-transcriptomics-heavy-processing
```

## Step 4: Create IAM Roles

Create IAM roles for Lambda and Batch:

```bash
# Create Lambda execution role
aws iam create-role \
    --role-name spatial-transcriptomics-lambda-role \
    --assume-role-policy-document file://infrastructure/policies/lambda-trust-policy.json

# Attach policies to Lambda role
aws iam attach-role-policy \
    --role-name spatial-transcriptomics-lambda-role \
    --policy-arn arn:aws:iam::aws:policy/service-role/AWSLambdaBasicExecutionRole

# Create custom policy for Lambda
aws iam create-policy \
    --policy-name spatial-transcriptomics-lambda-policy \
    --policy-document file://infrastructure/policies/lambda-policy.json

# Attach custom policy to Lambda role
aws iam attach-role-policy \
    --role-name spatial-transcriptomics-lambda-role \
    --policy-arn arn:aws:iam::account-id:policy/spatial-transcriptomics-lambda-policy

# Create Batch job role
aws iam create-role \
    --role-name spatial-transcriptomics-batch-job-role \
    --assume-role-policy-document file://infrastructure/policies/batch-trust-policy.json

# Create custom policy for Batch
aws iam create-policy \
    --policy-name spatial-transcriptomics-batch-policy \
    --policy-document file://infrastructure/policies/batch-policy.json

# Attach custom policy to Batch role
aws iam attach-role-policy \
    --role-name spatial-transcriptomics-batch-job-role \
    --policy-arn arn:aws:iam::account-id:policy/spatial-transcriptomics-batch-policy
```

## Step 5: Create Lambda Functions

Package and deploy the Lambda functions:

```bash
# Package job dispatcher
cd backend/lambda_functions/job_dispatcher
pip install -r requirements.txt -t .
zip -r function.zip .

# Deploy job dispatcher
aws lambda create-function \
    --function-name spatial-transcriptomics-job-dispatcher \
    --runtime python3.9 \
    --role arn:aws:iam::account-id:role/spatial-transcriptomics-lambda-role \
    --handler function.lambda_handler \
    --zip-file fileb://function.zip \
    --environment Variables="{DYNAMODB_TABLE=spatial-transcriptomics-jobs,BATCH_JOB_QUEUE=spatial-transcriptomics-job-queue,BATCH_JOB_DEFINITION_PREFIX=spatial-transcriptomics-}"

# Package job status updater
cd ../job_status_updater
pip install -r requirements.txt -t .
zip -r function.zip .

# Deploy job status updater
aws lambda create-function \
    --function-name spatial-transcriptomics-job-status-updater \
    --runtime python3.9 \
    --role arn:aws:iam::account-id:role/spatial-transcriptomics-lambda-role \
    --handler function.lambda_handler \
    --zip-file fileb://function.zip \
    --environment Variables="{DYNAMODB_TABLE=spatial-transcriptomics-jobs}"
```

Configure event sources for the Lambda functions:

```bash
# Configure SQS trigger for job dispatcher
aws lambda create-event-source-mapping \
    --function-name spatial-transcriptomics-job-dispatcher \
    --event-source-arn arn:aws:sqs:region:account-id:spatial-transcriptomics-light-processing \
    --batch-size 1

aws lambda create-event-source-mapping \
    --function-name spatial-transcriptomics-job-dispatcher \
    --event-source-arn arn:aws:sqs:region:account-id:spatial-transcriptomics-heavy-processing \
    --batch-size 1

# Create CloudWatch Events rule for job status updater
aws events put-rule \
    --name spatial-transcriptomics-status-checker \
    --schedule-expression "rate(2 minutes)"

# Add permission for CloudWatch Events
aws lambda add-permission \
    --function-name spatial-transcriptomics-job-status-updater \
    --statement-id CloudWatchEventsPermission \
    --action lambda:InvokeFunction \
    --principal events.amazonaws.com \
    --source-arn arn:aws:events:region:account-id:rule/spatial-transcriptomics-status-checker

# Configure CloudWatch Events target
aws events put-targets \
    --rule spatial-transcriptomics-status-checker \
    --targets Id=1,Arn=arn:aws:lambda:region:account-id:function:spatial-transcriptomics-job-status-updater
```

## Step 6: Set Up AWS Batch

Create a Compute Environment:

```bash
aws batch create-compute-environment \
    --compute-environment-name spatial-transcriptomics-compute-env \
    --type MANAGED \
    --state ENABLED \
    --compute-resources type=SPOT,minvCpus=0,maxvCpus=64,desiredvCpus=0,instanceTypes=c5.large,c5.xlarge,c5.2xlarge,r5.large,r5.xlarge,subnets=subnet-xxxx,securityGroupIds=sg-xxxx,instanceRole=ecsInstanceRole \
    --service-role AWSBatchServiceRole
```

Create a Job Queue:

```bash
aws batch create-job-queue \
    --job-queue-name spatial-transcriptomics-job-queue \
    --state ENABLED \
    --priority 1 \
    --compute-environment-order order=1,computeEnvironment=spatial-transcriptomics-compute-env
```

## Step 7: Build and Push Docker Images

Build the base Docker image:

```bash
cd infrastructure/docker
docker build -t spatial-transcriptomics-base:latest -f base.Dockerfile .
```

Create ECR repositories:

```bash
aws ecr create-repository --repository-name spatial-transcriptomics-normalization
aws ecr create-repository --repository-name spatial-transcriptomics-dimensionality-reduction
aws ecr create-repository --repository-name spatial-transcriptomics-clustering
aws ecr create-repository --repository-name spatial-transcriptomics-differential-expression
```

Build and push job-specific Docker images:

```bash
# Get ECR login
aws ecr get-login-password | docker login --username AWS --password-stdin account-id.dkr.ecr.region.amazonaws.com

# Build and push normalization image
cd backend/batch_jobs/normalization
docker build -t account-id.dkr.ecr.region.amazonaws.com/spatial-transcriptomics-normalization:latest .
docker push account-id.dkr.ecr.region.amazonaws.com/spatial-transcriptomics-normalization:latest

# Repeat for other job types
```

## Step 8: Create Batch Job Definitions

Create job definitions for each analysis type:

```bash
aws batch register-job-definition \
    --job-definition-name spatial-transcriptomics-normalization \
    --type container \
    --container-properties '{"image":"account-id.dkr.ecr.region.amazonaws.com/spatial-transcriptomics-normalization:latest","vcpus":2,"memory":8000,"jobRoleArn":"arn:aws:iam::account-id:role/spatial-transcriptomics-batch-job-role","environment":[{"name":"DYNAMODB_TABLE","value":"spatial-transcriptomics-jobs"}]}'

aws batch register-job-definition \
    --job-definition-name spatial-transcriptomics-dimensionality_reduction \
    --type container \
    --container-properties '{"image":"account-id.dkr.ecr.region.amazonaws.com/spatial-transcriptomics-dimensionality-reduction:latest","vcpus":4,"memory":16000,"jobRoleArn":"arn:aws:iam::account-id:role/spatial-transcriptomics-batch-job-role","environment":[{"name":"DYNAMODB_TABLE","value":"spatial-transcriptomics-jobs"}]}'

aws batch register-job-definition \
    --job-definition-name spatial-transcriptomics-clustering \
    --type container \
    --container-properties '{"image":"account-id.dkr.ecr.region.amazonaws.com/spatial-transcriptomics-clustering:latest","vcpus":4,"memory":16000,"jobRoleArn":"arn:aws:iam::account-id:role/spatial-transcriptomics-batch-job-role","environment":[{"name":"DYNAMODB_TABLE","value":"spatial-transcriptomics-jobs"}]}'

aws batch register-job-definition \
    --job-definition-name spatial-transcriptomics-differential_expression \
    --type container \
    --container-properties '{"image":"account-id.dkr.ecr.region.amazonaws.com/spatial-transcriptomics-differential-expression:latest","vcpus":8,"memory":32000,"jobRoleArn":"arn:aws:iam::account-id:role/spatial-transcriptomics-batch-job-role","environment":[{"name":"DYNAMODB_TABLE","value":"spatial-transcriptomics-jobs"}]}'
```

## Step 9: Update Environment Variables

Update your `.env` file with the job definition ARNs:

```
BATCH_NORMALIZATION_JOB_DEF=spatial-transcriptomics-normalization
BATCH_DIM_REDUCTION_JOB_DEF=spatial-transcriptomics-dimensionality_reduction
BATCH_CLUSTERING_JOB_DEF=spatial-transcriptomics-clustering
BATCH_DIFF_EXP_JOB_DEF=spatial-transcriptomics-differential_expression
DYNAMODB_JOB_TABLE=spatial-transcriptomics-jobs
```

## Step 10: Test the Infrastructure

Run a simple test to verify the infrastructure is working:

```bash
# Create a test message
echo '{"job_id":"test-job","job_type":"normalization","input_data":{"s3_uri":"s3://spatial-transcriptomics-data/test/input.h5ad","filename":"test"},"job_params":{"normalize":true},"submitted_at":"2023-01-01T00:00:00Z","status":"SUBMITTED"}' > test_message.json

# Send message to SQS
aws sqs send-message \
    --queue-url https://sqs.region.amazonaws.com/account-id/spatial-transcriptomics-light-processing \
    --message-body file://test_message.json

# Check DynamoDB for job status
aws dynamodb get-item \
    --table-name spatial-transcriptomics-jobs \
    --key '{"job_id":{"S":"test-job"}}'
```

## Cleanup Instructions

If you need to remove all created resources:

```bash
# Delete Lambda functions
aws lambda delete-function --function-name spatial-transcriptomics-job-dispatcher
aws lambda delete-function --function-name spatial-transcriptomics-job-status-updater

# Delete SQS queues
aws sqs delete-queue --queue-url https://sqs.region.amazonaws.com/account-id/spatial-transcriptomics-light-processing
aws sqs delete-queue --queue-url https://sqs.region.amazonaws.com/account-id/spatial-transcriptomics-heavy-processing

# Delete Batch job definitions
aws batch deregister-job-definition --job-definition spatial-transcriptomics-normalization:1
aws batch deregister-job-definition --job-definition spatial-transcriptomics-dimensionality_reduction:1
aws batch deregister-job-definition --job-definition spatial-transcriptomics-clustering:1
aws batch deregister-job-definition --job-definition spatial-transcriptomics-differential_expression:1

# Delete Batch job queue
aws batch update-job-queue --job-queue spatial-transcriptomics-job-queue --state DISABLED
aws batch delete-job-queue --job-queue spatial-transcriptomics-job-queue

# Delete Batch compute environment
aws batch update-compute-environment --compute-environment spatial-transcriptomics-compute-env --state DISABLED
aws batch delete-compute-environment --compute-environment spatial-transcriptomics-compute-env

# Delete DynamoDB table
aws dynamodb delete-table --table-name spatial-transcriptomics-jobs

# Delete S3 bucket (empty it first)
aws s3 rm s3://spatial-transcriptomics-data --recursive
aws s3api delete-bucket --bucket spatial-transcriptomics-data

# Delete IAM roles and policies
aws iam detach-role-policy --role-name spatial-transcriptomics-lambda-role --policy-arn arn:aws:iam::aws:policy/service-role/AWSLambdaBasicExecutionRole
aws iam detach-role-policy --role-name spatial-transcriptomics-lambda-role --policy-arn arn:aws:iam::account-id:policy/spatial-transcriptomics-lambda-policy
aws iam delete-policy --policy-arn arn:aws:iam::account-id:policy/spatial-transcriptomics-lambda-policy
aws iam delete-role --role-name spatial-transcriptomics-lambda-role

aws iam detach-role-policy --role-name spatial-transcriptomics-batch-job-role --policy-arn arn:aws:iam::account-id:policy/spatial-transcriptomics-batch-policy
aws iam delete-policy --policy-arn arn:aws:iam::account-id:policy/spatial-transcriptomics-batch-policy
aws iam delete-role --role-name spatial-transcriptomics-batch-job-role
```