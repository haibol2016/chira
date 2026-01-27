# ChiRA Dockerfile
# Uses micromamba for lightweight conda package management
# Designed for compatibility with both Docker and Singularity/Apptainer

FROM mambaorg/micromamba:1.5.8-jammy

# Create non-root user for Singularity compatibility
# Use high-range UID to avoid conflicts with host users
USER root
RUN groupadd -r -g 9001 chira && \
    useradd -r -u 9001 -g chira -m -d /home/chira chira && \
    mkdir -p /app /home/chira/{data,output,scratch} && \
    chown -R chira:chira /app /home/chira

WORKDIR /app

# Install system dependencies and set up conda in single layer
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        curl \
        wget \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    # Install conda packages
    micromamba install -y -n base -c conda-forge -c bioconda \
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
        && \
    micromamba clean -afy && \
    # Fix permissions for conda directory
    chown -R chira:chira /opt/conda

# Set environment variables for both Docker and Singularity
ENV CONDA_PREFIX=/opt/conda
ENV PATH=/app:/opt/conda/bin:${PATH}
ENV PYTHONPATH=/app
ENV HOME=/home/chira
ENV TMPDIR=/tmp

# Create Singularity environment script (automatically sourced by Singularity)
RUN mkdir -p /.singularity.d/env && \
    echo '#!/bin/bash' > /.singularity.d/env/91-environment.sh && \
    echo 'export CONDA_PREFIX=/opt/conda' >> /.singularity.d/env/91-environment.sh && \
    echo 'export PATH=/app:/opt/conda/bin:${PATH}' >> /.singularity.d/env/91-environment.sh && \
    echo 'export PYTHONPATH=/app' >> /.singularity.d/env/91-environment.sh && \
    echo 'export HOME=/home/chira' >> /.singularity.d/env/91-environment.sh && \
    echo 'export TMPDIR=/tmp' >> /.singularity.d/env/91-environment.sh && \
    chmod +x /.singularity.d/env/91-environment.sh

# Copy the ChiRA codebase and set permissions
COPY --chown=chira:chira *.py LICENSE ./

# Make Python scripts executable
RUN chmod +x *.py && \
    chown -R chira:chira /app

# Switch to non-root user
USER chira


# Default command
CMD ["/bin/bash"]

