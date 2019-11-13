FROM gliderlabs/alpine:3.9

RUN apk --no-cache --update-cache -U add tzdata && \
    cp /usr/share/zoneinfo/America/Chicago /etc/localtime && \
    apk --no-cache --update-cache add \
        bash \
        curl \
        openjdk8 && \
     rm -rf /var/cache/apk/*

WORKDIR /opt

RUN curl -L https://downloads.sourceforge.net/project/bbmap/BBMap_38.72.tar.gz -o BBMap.tar.gz && \
    tar xzf BBMap.tar.gz && \
    mv bbmap /usr/local && \
    rm -rf BBMap.tar.gz

ENV PATH="/usr/local/bbmap:${PATH}"

CMD ["bash", "bbduk.sh"]