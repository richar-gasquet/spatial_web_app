FROM spatial-transcriptomics-base:latest

# Copy job-specific scripts
COPY process.py /app/

# Set entry point
ENTRYPOINT ["python", "/app/process.py"]