FROM mambaorg/micromamba:git-86a1dd8-cuda13.2.1-ubuntu24.04

ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root

# Install system dependencies
RUN micromamba install -y -c bioconda -c conda-forge \
    samtools \
    picard \
    perl \
    python=3.10

# Set pip mirror
RUN pip config set global.index-url https://mirrors.cloud.tencent.com/pypi/simple

# Copy and install pindel-tool package
COPY . /app
WORKDIR /app
RUN pip install .

# Clean up
RUN micromamba clean --all --yes

# # Set environment
# ENV PATH=/opt/conda/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/conda/lib:$LD_LIBRARY_PATH

# Set entrypoint
# ENTRYPOINT ["/opt/conda/bin/pindel-tool"]
# CMD ["--help"]
