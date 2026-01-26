# ChiRA Dockerfile
# Uses micromamba for lightweight conda package management

FROM mambaorg/micromamba:1.5.8-jammy

# Set working directory
WORKDIR /app

# Switch to root user for system package installation
USER root

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Switch to micromamba user
USER $MAMBA_USER

# Set up conda channels and install packages
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.11 \
    biopython \
    bcbiogff \
    pysam \
    requests \
    pyliftover \
    bwa \
    samtools \
    bedtools \
    intarna \
    && micromamba clean -afy

# Note: CLAN and blockbuster are not available in conda-forge/bioconda
# They need to be installed separately if required

# Activate the base environment
ENV CONDA_PREFIX=/opt/conda
ENV PATH=/app:$CONDA_PREFIX/bin:$PATH

# Copy the ChiRA codebase
COPY --chown=$MAMBA_USER:$MAMBA_USER *.py ./
#COPY --chown=$MAMBA_USER:$MAMBA_USER *.md ./
COPY --chown=$MAMBA_USER:$MAMBA_USER LICENSE ./

# Make Python scripts executable
# List files explicitly to avoid glob expansion issues
USER root
RUN chmod +x *.py
USER $MAMBA_USER

# Set Python path to include current directory
ENV PYTHONPATH=/app

# Default command
CMD ["/bin/bash"]

