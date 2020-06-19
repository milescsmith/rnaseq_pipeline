FROM continuumio/miniconda3:latest

LABEL maintainer="miles-smith@omrf.org"

ENV APT_LISTCHANGES_FRONTEND=none
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update && \
   apt-get dist-upgrade -u -y && \
   apt-get -y install ttf-dejavu && \
   apt-get clean && \
   rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

COPY conda.yml .
RUN \
   conda update conda && \
   conda env update -n root -f conda.yml && \
   conda clean -a
