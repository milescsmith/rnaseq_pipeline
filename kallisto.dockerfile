FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
		build-essential \
		cmake \
		python \
		python-pip \
		python-dev \
		hdf5-tools \
		libhdf5-dev \
		hdf5-helpers \
		libhdf5-serial-dev \
		git \
		apt-utils

WORKDIR /root
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.45.0/kallisto_linux-v0.45.0.tar.gz && \
    tar xzf kallisto_linux-v0.45.0.tar.gz && \
    mv kallisto_linux-v0.45.0/kallisto /usr/bin/ && \
    rm -rf kallisto_linux-v0.45.0 kallisto_linux-v0.45.0.tar.gz

CMD /usr/bin/kallisto