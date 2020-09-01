FROM gliderlabs/alpine:3.9

RUN apk --no-cache --update-cache -U add tzdata && \
    cp /usr/share/zoneinfo/America/Chicago /etc/localtime && \
    apk --no-cache --update-cache add \
        curl \
        openjdk8 \
        perl && \
     rm -rf /var/cache/apk/*

WORKDIR /opt

RUN curl -L https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip -o fastqc.zip && \
    unzip fastqc.zip && \
    chmod 750 FastQC/fastqc && \
    ln -s FastQC/fastqc /usr/local/bin/fastqc

WORKDIR /root

ENV PATH="/opt/FastQC:${PATH}"

CMD ["fastqc", "--help"]