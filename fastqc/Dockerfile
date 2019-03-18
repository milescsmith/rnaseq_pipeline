FROM ubuntu:18.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        wget \
        hdf5-tools \
        libhdf5-dev \
        hdf5-helpers \
        libhdf5-serial-dev \
        git \
        apt-utils \
        python3 \
        python3-pip \
        unzip \
        openjdk-11-jre \
        fastqc && \
    apt-get clean && \
    rm -rf /tmp/downloaded_packages/ /tmp/*.rds && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install setuptools wheel && pip3 install kubernetes

CMD ["fastqc"]
