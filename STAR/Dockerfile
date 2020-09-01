FROM gliderlabs/alpine:3.7

LABEL author="Miles Smith miles-smith@omrf.org"

RUN apk --no-cache --update-cache -U add tzdata && \
    cp /usr/share/zoneinfo/America/Chicago /etc/localtime && \
    apk --no-cache --update-cache add \
        curl

WORKDIR /opt

RUN \
    curl -L https://github.com/alexdobin/STAR/archive/2.7.5c.tar.gz -o 2.7.5c.tar.gz && \
    tar -xzf 2.7.5c.tar.gz && \
    mv STAR-2.7.5c/bin/Linux_x86_64_static/* /usr/local/bin && \
    rm -fr 2.7.5c.tar.gz && \
    rm -fr STAR-2.7.5c

CMD ["STAR", "--help"]