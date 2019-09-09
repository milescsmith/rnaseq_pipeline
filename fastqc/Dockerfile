FROM frolvlad/alpine-java

RUN apk --no-cache --update-cache add \
    curl \
    perl && \
    rm -rf /var/cache/apk/*
    
WORKDIR /root
RUN curl -o fastqc.zip -L https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
    unzip fastqc.zip && \
    mv FastQC /usr/local/FastQC && \
    rm -fr fastqc.zip && \
    chmod 755 /usr/local/FastQC/fastqc && \
    ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc

CMD ["-help"]
ENTRYPOINT ["/usr/local/bin/fastqc"]