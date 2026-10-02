
FROM continuumio/miniconda3:24.7.1-0
 
# Avoids apt hanging on interactive prompts mid-build.
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /bulk_rna_seq
 
# gsutil — needed if any task in this image reads/writes GCS paths directly.
RUN apt-get update && apt-get install -y \
        curl \
        git \
        gnupg \
        build-essential \
        python3-dev && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" \
        | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get update && apt-get install -y google-cloud-cli && \
    rm -rf /var/lib/apt/lists/*
 

RUN git clone --depth 1 https://github.com/klarman-cell-observatory/bulk_rna_workflow.git /bulk_rna_seq/repo
 
COPY bulkrna_env.yml /bulk_rna_seq/bulkrna_env.yml
 
RUN conda install -n base -c conda-forge mamba -y && \
    mamba env create -f /bulk_rna_seq/bulkrna_env.yml -n bulk_rna && \
    conda clean -afy
 
ENV PATH=/opt/conda/envs/bulk_rna/bin:$PATH
 
RUN Rscript -e 'library(DESeq2); library(apeglm); library(ggplot2); library(pheatmap); library(RColorBrewer); library(ggrepel); library(stringr); cat("R OK\n")'
RUN python -c "import pandas, numpy, matplotlib, openpyxl, adjustText, firecloud.api; print('python OK')"
RUN gsutil --version
RUN ls -la /bulk_rna_seq/repo/scripts/