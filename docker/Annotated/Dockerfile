FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

LABEL maintainer="Servicios Inmegen"
LABEL version="Annotated 1.0"

# Be sure to symlink /var/lib/docker to a larger disk.
# Build the image: docker build -t pipelinesinmegen/pipelines_inmegen:snpeff -f Dockerfile .

# Set shell
RUN set -e && ln -sf bash /bin/sh

### Update and install linux tools 
RUN apt-get update && \
    apt-get install -y \
    ant \
    automake \
    apt-utils \
    build-essential \
    bzip2 \
    ca-certificates \
    curl \
    gawk \
    gcc \
    git \
    g++ \
    hmmer \
    less \
    locales \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libpng-dev \
    libpcre2-dev \
    libpcre3 \
    libpcre3-dev \
    libncurses5-dev \
    libnss-sss \
    libcairo2-dev \
    libxt-dev \
    libfontconfig1-dev \
    libudunits2-dev \
    libgdal-dev \
    less \
    locales \
    make \
    nano \
    nim \
    openjdk-21-jdk \
    parallel \
    perl \
    pigz \
    pkg-config \
    python3-dev \
    python3-numpy \
    python3-setuptools \
    python3-pip \
    software-properties-common \
    tabix \
    tar \
    tzdata \
    unzip \
    wget \
    zip \
    zlib1g \
    zlib1g-dev && \
    apt-get -y autoremove &&\
    apt-get clean &&\
    rm -rf /var/lib/apt/lists/*

#### Download and install Annovar
RUN wget Annovar-URL && \
    tar -zxvf annovar.latest.tar.gz && \
    rm annovar.latest.tar.gz

### Wrappers for Annovar annotate_variation.pl
RUN mkdir -p /annovar/bin && \
    echo '#!/usr/bin/env bash' > /annovar/bin/table_annovar && \
    echo 'perl /annovar/table_annovar.pl ${@}' >> /annovar/bin/table_annovar && \
    echo '#!/usr/bin/env bash' > /annovar/bin/convert2annovar && \
    echo 'perl /annovar/convert2annovar.pl ${@}' >> /annovar/bin/convert2annovar && \
    echo '#!/usr/bin/env bash' > /annovar/bin/annotate_variation && \
    echo 'perl /annovar/annotate_variation.pl ${@}' >> /annovar/bin/annotate_variation && \
    chmod a+x /annovar/bin/*

### Download and install SnpEff
RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip && \
    chmod a+x /snpEff/snpEff.jar && \
    java -jar /snpEff/snpEff.jar download GRCh38.99 && \
    rm snpEff_latest_core.zip

### Wrappers for snpEff, SnpSift and ClinEff 
RUN mkdir -p /snpEff/bin && \
    echo '#!/usr/bin/env bash' > /snpEff/bin/snpEff && \
    echo 'java -jar /snpEff/snpEff.jar ${@}' >> /snpEff/bin/snpEff && \
    echo '#!/usr/bin/env bash' > /snpEff/bin/SnpSift && \
    echo 'java -jar /snpEff/SnpSift.jar ${@}' >> /snpEff/bin/SnpSift && \
    chmod a+x /snpEff/bin/*

### Environment variables
ENV TZ "America/Mexico_City"
ENV PATH /annovar/bin/:/snpEff/bin:${PATH}

### Command run
CMD ["/bin/bash"]
