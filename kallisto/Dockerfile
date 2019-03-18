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
        python3-pip && \
    apt-get clean && \
    rm -rf /tmp/downloaded_packages/ /tmp/*.rds && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /root
RUN wget --no-check-certificate https://github.com/pachterlab/kallisto/releases/download/v0.45.0/kallisto_linux-v0.45.0.tar.gz && \
    tar xzf kallisto_linux-v0.45.0.tar.gz && \
    mv kallisto_linux-v0.45.0/kallisto /usr/bin/ && \
    rm -rf kallisto_linux-v0.45.0 kallisto_linux-v0.45.0.tar.gz

RUN pip3 install setuptools wheel && pip3 install kubernetes

CMD /usr/bin/kallisto
