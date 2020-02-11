FROM frolvlad/alpine-gcc

LABEL version="1.0"
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
        linux-headers \
        && \
    rm -rf /tmp/* /var/cache/apk/*

RUN pip3 install --no-cache-dir snakemake

CMD ["snakemake", "--help"]