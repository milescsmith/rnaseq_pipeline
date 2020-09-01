FROM quay.io/snakemake/snakemake:latest

LABEL description="Tools to run an RNAseq QC and alignemtns processing Snakemake pipeline in Google Kubenetes Engine"
LABEL tags="Genomics"
LABEL author="Miles Smith <miles-smith@omrf.org>"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    hdf5-tools \
    libhdf5-dev \
    hdf5-helpers \
    libhdf5-serial-dev \
    git \
    apt-utils \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

WORKDIR /root

RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda install -y -c conda-forge \
    hdf5 \
    kubernetes=1.13.4 \
    && conda install -y -c bioconda \
    snakemake=5.4.3 \
    multiqc=1.7 \
    kallisto=0.45.0 \
    bbmap=38.22 \
    fastqc=0.11.8 \
    star=2.7.0e \
    samtools=1.9 \
    stringtie=1.3.5 \
    && conda clean -y --all
