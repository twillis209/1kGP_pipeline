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

# Base conda env: snakemake + python
RUN mamba install -y \
        -c conda-forge \
        -c bioconda \
        "snakemake-minimal>=9.6" \
        "python>=3.11,<3.14" \
    && mamba clean -afy

# Pre-bake all conda environments
RUN mkdir -p /tmp/envs && \
    printf 'channels:\n  - conda-forge\n  - bioconda\ndependencies:\n  - pandas\n  - polars\n' \
        > /tmp/envs/global.yaml && \
    printf 'channels:\n  - conda-forge\n  - bioconda\ndependencies:\n  - r-base\n  - r-data.table >= 1.17\n  - r-ggplot2\n' \
        > /tmp/envs/r.yaml && \
    printf 'name: ldsc\nchannels:\n  - conda-forge\n  - bioconda\ndependencies:\n  - python>=3.9,<3.14\n  - numpy\n  - scipy\n  - pandas\n  - bitarray\n  - pybedtools\n  - setuptools\n  - pip\n' \
        > /tmp/envs/ldsc.yaml && \
    mamba env create -p /opt/snakemake-envs/global --file /tmp/envs/global.yaml && \
    mamba env create -p /opt/snakemake-envs/r     --file /tmp/envs/r.yaml     && \
    mamba env create -p /opt/snakemake-envs/ldsc  --file /tmp/envs/ldsc.yaml  && \
    /opt/snakemake-envs/ldsc/bin/pip install --no-deps \
        git+https://github.com/bulik/ldsc.git && \
    mamba clean -afy && \
    rm -rf /tmp/envs

# plink2 and plink binaries
RUN wget -q -O /tmp/plink2.zip \
        https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20260310.zip \
    && unzip -j /tmp/plink2.zip plink2 -d /usr/local/bin \
    && chmod +x /usr/local/bin/plink2 \
    && rm /tmp/plink2.zip \
    && wget -q -O /tmp/plink.zip \
        https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip \
    && unzip -j /tmp/plink.zip plink -d /usr/local/bin \
    && chmod +x /usr/local/bin/plink \
    && rm /tmp/plink.zip

ENV PATH="/opt/conda/bin:$PATH"
ENV CONDA_DEFAULT_ENV=base
