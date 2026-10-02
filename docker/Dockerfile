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

# Pull the repo's scripts/ folder into the image at a fixed path, so the analysis scripts are
# available in-container even on a run that doesn't pass build_count_matrix_script/qc_plots_script/
# run_deseq2_script as WDL File inputs.
# TODO: confirmed klarman-cell-observatory/bulk_rna_workflow exists and has a scripts/ folder
# (verified earlier in this conversation); your original Dockerfile pointed at bulk_rna_seq, which
# I couldn't find under that org — swap the URL below if that name is actually correct.
RUN git clone --depth 1 https://github.com/klarman-cell-observatory/bulk_rna_workflow.git /bulk_rna_seq/repo

# bulkrna_env.yml is your real env from the notebook (renamed to bulk_rna here, was amyz_bulkrna_env), plus DESeq2/apeglm/
# pheatmap/ggrepel/pandas/openpyxl/adjusttext added for the three R scripts + the WDL's python
# glue — see the "added" block near the top of that file for exactly what's new vs. what you
# already had (fastqc/multiqc/tidyverse/jupyter kernels etc. are all untouched).
COPY bulkrna_env.yml /bulk_rna_seq/bulkrna_env.yml

RUN conda env create -f /bulk_rna_seq/bulkrna_env.yml -n bulk_rna && \
    conda clean -afy

# Putting the env's bin/ first on PATH is enough for Rscript/python3 to resolve correctly inside
# a WDL task's command block — no `conda activate` needed (it wouldn't persist into a fresh task
# shell anyway).
ENV PATH=/opt/conda/envs/bulk_rna/bin:$PATH

# Fail the build here, at build time, rather than three hours into a Terra submission.
RUN Rscript -e 'library(DESeq2); library(apeglm); library(ggplot2); library(pheatmap); library(RColorBrewer); library(ggrepel); library(stringr); cat("R OK\n")'
RUN python -c "import pandas, numpy, matplotlib, openpyxl, adjustText, firecloud.api; print('python OK')"
RUN gsutil --version
RUN ls -la /bulk_rna_seq/repo/scripts/