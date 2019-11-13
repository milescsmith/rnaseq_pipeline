FROM gliderlabs/alpine:3.9

RUN apk --no-cache --update-cache -U add tzdata && \
    cp /usr/share/zoneinfo/America/Chicago /etc/localtime && \
    apk --no-cache --update-cache add \
        curl && \
     rm -rf /var/cache/apk/*

WORKDIR /opt

RUN curl -k -L https://downloads.sourceforge.net/project/subread/subread-2.0.0/subread-2.0.0-Linux-x86_64.tar.gz -o subread.tar.gz && \
    tar xzf subread.tar.gz && \
    rm -rf subread.tar.gz

ENV PATH="/opt/subread-2.0.0-Linux-x86_64/bin:${PATH}"

CMD ["featureCounts", "--help"]