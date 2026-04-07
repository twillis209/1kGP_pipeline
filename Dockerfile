# syntax=docker/dockerfile:1
FROM condaforge/mambaforge:24.3.0-0

LABEL org.opencontainers.image.title="1kGP pipeline"
LABEL org.opencontainers.image.description="Snakemake pipeline for 1000 Genomes Project data processing"

# System packages
RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        unzip \
        zstd \
        git \
        ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Base conda env: python
RUN mamba install -y \
        -c conda-forge \
        "python>=3.11,<3.14" \
    && mamba clean -afy

# Pre-bake all conda environments
COPY workflow/envs/global.yaml /tmp/envs/global.yaml
COPY workflow/envs/r.yaml      /tmp/envs/r.yaml
COPY workflow/envs/ldsc.yaml   /tmp/envs/ldsc.yaml

RUN mamba env create -p /opt/snakemake-envs/global --file /tmp/envs/global.yaml && \
    mamba env create -p /opt/snakemake-envs/r     --file /tmp/envs/r.yaml     && \
    mamba env create -p /opt/snakemake-envs/ldsc  --file /tmp/envs/ldsc.yaml  && \
    /opt/snakemake-envs/ldsc/bin/pip install --no-deps \
        git+https://github.com/bulik/ldsc.git && \
    mamba clean -afy && \
    rm -rf /tmp/envs

# plink2 and plink binaries
RUN wget -q -O /tmp/plink2.zip \
        https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20260311.zip \
    && unzip -j /tmp/plink2.zip plink2 -d /usr/local/bin \
    && chmod +x /usr/local/bin/plink2 \
    && rm /tmp/plink2.zip \
    && wget -q -O /tmp/plink.zip \
        https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip \
    && unzip -j /tmp/plink.zip plink -d /usr/local/bin \
    && chmod +x /usr/local/bin/plink \
    && rm /tmp/plink.zip

ENV CONDA_DEFAULT_ENV=base
