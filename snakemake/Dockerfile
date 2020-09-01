FROM frolvlad/alpine-gcc

LABEL version="1.1"
LABEL author="Miles Smith <miles-smith@omrf.org>"
LABEL maintainer="Miles Smith <miles-smith@omrf.org>"
LABEL description="Container for running dockerized Snakemake scripts"

RUN apk add --no-cache --virtual=.build-dependencies \
        bash \
        grep \
        sed \
        curl \
        docker \
        python3 \
        python3-dev \
        py3-pip \
        linux-headers \
        && \
    rm -rf /tmp/* /var/cache/apk/*

RUN pip3 install --no-cache-dir snakemake==5.23.0

CMD ["snakemake", "--help"]