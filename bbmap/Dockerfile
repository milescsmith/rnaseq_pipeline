FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
    build-essential \
    openjdk-11-jre \
    wget \
    && apt-get clean \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /root
RUN wget -O BBMap.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_38.56.tar.gz/download && \
    tar xzf BBMap.tar.gz && \
    mv bbmap /usr/local && \
    rm -rf BBMap.tar.gz

ENV PATH="/usr/local/bbmap:${PATH}"

CMD ["bbduk.sh"]
