FROM continuumio/miniconda3:latest

WORKDIR /bulk_rna_seq

RUN mkdir -p /bulk_rna_seq/data
RUN mkdir -p /bulk_rna_seq/results

RUN git clone https://github.com/klarman-cell-observatory/bulk_rna_seq /bulk_rna_seq_main

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

RUN conda env create -f bulkrna_env.yml -n bulk_rna

RUN CONDA_DEFAULT_ENV=bulk_rna
ENV PATH=/opt/conda/envs/bulk_rna/bin:$PATH

RUN /bin/bash -c "source activate bulk_rna"