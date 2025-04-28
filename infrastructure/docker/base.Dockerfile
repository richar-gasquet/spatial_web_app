FROM python:3.11-bullseye

# Install system dependencies
RUN apt-get update && apt-get install -y  \
    gcc \
    g++ \
    libhdf5-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Set up working directory
WORKDIR /app

# Install Python dependencies for spatial transcriptomics analysis
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy common scripts
COPY common/ /app/common/

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PYTHONPATH=/app

# Entry point will be specified in derived images
CMD ["echo", "Base image - use a specific job image instead"]