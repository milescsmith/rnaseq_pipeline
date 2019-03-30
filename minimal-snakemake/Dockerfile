FROM python:3.6-alpine

RUN apk add --no-cache --virtual=.build-dependencies wget ca-certificates bash python3 gcc musl-dev linux-headers libffi-dev libressl-dev \
    && rm -rf /tmp/* /var/cache/apk/*
RUN pip install snakemake kubernetes google-cloud-sdk
CMD ["snakemake"]
