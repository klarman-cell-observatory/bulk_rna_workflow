FROM continuumio/miniconda3:24.7.1-0

# Avoids apt hanging on interactive prompts mid-build.
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /bulk_rna_seq

# gsutil — needed if any task in this image reads/writes GCS paths directly.
RUN apt-get update && apt-get install -y \
        curl \
        gnupg \
        build-essential \
        python3-dev && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" \
        | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get update && apt-get install -y google-cloud-sdk && \
    rm -rf /var/lib/apt/lists/*

# bulkrna_env.yml is your real env (amyz_bulkrna_env, from the notebook), plus DESeq2/apeglm/
# pheatmap/ggrepel/pandas/openpyxl/adjusttext added for the three R scripts + the WDL's python
# glue — see the "added" block near the top of that file for exactly what's new vs. what you
# already had (fastqc/multiqc/tidyverse/jupyter kernels etc. are all untouched).
COPY bulkrna_env.yml /bulk_rna_seq/bulkrna_env.yml

# ~300 packages, most pinned to exact versions, now with a few new unpinned additions layered in
# — that's a big solve for classic conda. mamba's solver handles this in a fraction of the time
# and is far less likely to hang or fail outright.
RUN conda install -n base -c conda-forge mamba -y && \
    mamba env create -f /bulk_rna_seq/bulkrna_env.yml -n amyz_bulkrna_env && \
    conda clean -afy

# Putting the env's bin/ first on PATH is enough for Rscript/python3 to resolve correctly inside
# a WDL task's command block — no `conda activate` needed (it wouldn't persist into a fresh task
# shell anyway).
ENV PATH=/opt/conda/envs/amyz_bulkrna_env/bin:$PATH

# Fail the build here, at build time, rather than three hours into a Terra submission.
RUN Rscript -e 'library(DESeq2); library(apeglm); library(ggplot2); library(pheatmap); library(RColorBrewer); library(ggrepel); library(stringr); cat("R OK\n")'
RUN python -c "import pandas, numpy, matplotlib, openpyxl, adjustText; print('python OK')"
RUN gsutil --version