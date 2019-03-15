FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
	build-essential \
	openjdk-11-jre \
	wget \

ENV BBMAP_VERSION 38.20
ENV BBMAP_DIR /usr/local/bbmap

WORKDIR /root
RUN wget -O BBMap.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_38.42.tar.gz/download && \
    tar xzf BBMap_38.42.tar.gz && \
    mv bbmap /usr/local && \
    rm -rf BBMap_38.42.tar.gz

ENV PATH="/usr/local/bbmap:${PATH}"

CMD ["bbduk.sh"]