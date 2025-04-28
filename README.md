# Spatial Transcriptomics Analysis Pipeline

A multi-step processing application for spatial transcriptomics data with Streamlit frontend and AWS Batch backend for heavy computation.

## Overview

This application provides a user-friendly interface for analyzing spatial transcriptomics data through the following steps:

1. **Data Upload** - Upload spatial transcriptomics data files
2. **Data Preparation** - Calculate QC metrics and prepare data
3. **Data Visualization** - Visualize data distribution and QC metrics
4. **Quality Control & Filtering** - Filter spots/cells and genes
5. **Analysis** - Run various analyses including normalization, dimensionality reduction, clustering, and differential expression

The application uses AWS services to offload heavy computation to the cloud:
- AWS S3 for data storage
- AWS Batch for heavy computation
- AWS Lambda for job management
- DynamoDB for job tracking
- SQS for job queues

## Setup

1. Clone this repository
2. Create a virtual environment: `python -m venv venv`
3. Activate the virtual environment: `source venv/bin/activate` (Linux/Mac) or `venv\Scripts\activate` (Windows)
4. Install requirements: `pip install -r requirements.txt`
5. Copy `.env.example` to `.env` and fill in your AWS credentials and configuration
6. Deploy AWS resources (see AWS Setup below)
7. Run the application: `cd frontend && streamlit run app.py`

## AWS Setup

See the detailed AWS setup instructions in the documentation.

## Development

This project is structured as follows:

- `/frontend` - Streamlit application
- `/backend` - AWS Batch job scripts and Lambda functions
- `/docs` - Documentation

## License

[MIT License](LICENSE)