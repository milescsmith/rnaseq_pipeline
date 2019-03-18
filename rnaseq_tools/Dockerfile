FROM quay.io/snakemake/snakemake:latest

RUN apt-get update && apt-get install -y \
    build-essential \
    hdf5-tools \
    libhdf5-dev \
    hdf5-helpers \
    libhdf5-serial-dev \
    git \
    apt-utils

WORKDIR /root

run conda config --add channels conda-forge
run conda config --add channels bioconda
run conda install -c conda-forge hdf5
RUN conda install -c bioconda multiqc kallisto=0.45.0 bbmap=38.22 fastqc=0.11.8 star=2.7.0e samtools=1.9 stringtie=1.3.5 qualimap=2.2.2b picard=2.18.27 deeptools=3.2.0 qorts=1.3.0 
RUN conda install -c bioconda snakemake=5.4.3 
run conda install -c conda-forge kubernetes=1.13.4 google-cloud-sdk=238.0.0 google-cloud-storage=1.14.0 google-resumable-media=0.3.2