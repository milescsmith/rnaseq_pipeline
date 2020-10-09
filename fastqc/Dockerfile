FROM ubuntu:bionic

ENV DEBIAN_FRONTEND=noninteractive
RUN apt update --fix-missing && \
   apt install --no-install-recommends -y \
      curl \
      openjdk-11-jdk \
      openjdk-11-jre \
      perl \
      unzip && \
   apt-get clean && \
   rm -rf /tmp/downloaded_packages/* && \
   rm -rf /var/lib/apt/lists/*

WORKDIR /opt

RUN curl -L https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip -o fastqc.zip && \
    unzip fastqc.zip && \
    chmod 750 FastQC/fastqc && \
    ln -s FastQC/fastqc /usr/local/bin/fastqc

WORKDIR /root

ENV PATH="/opt/FastQC:${PATH}"

CMD ["fastqc", "--help"]