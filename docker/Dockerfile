FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

LABEL maintainer="Subdirección de genomica poblacional y subdirección de bioinformática del INMEGEN"
LABEL version="0.1"

# Building could takes several hours! Be sure to symlink /var/lib/docker to a larger disk.
# Build the image: docker build -t pipelines_inmegen -f Dockerfile .
# To run the container: docker run -i -t pipelines_inmegen /bin/bash

# Actualiza el sistema operativo e instala las herramientas adecuadas
RUN apt-get update && \
    apt-get install -y \
        ant \
        automake \
        apt-utils \
        build-essential \
        bzip2 \
        ca-certificates \
        curl \
        cutadapt \
        gawk \
        gcc \
        git \
        g++ \
        hmmer \
        kallisto \
        openjdk-8-jdk \
        openjdk-11-jdk \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libxml2-dev \
        libpng-dev \
        libpcre2-dev \
        libcairo2-dev \
        libxt-dev \
        libfontconfig1-dev \
        libudunits2-dev \
        libgdal-dev \
        less \
        locales \
        make \
        nano \
        parallel \
        perl \
        pigz \
        pkg-config \
        python3-dev \
        python3-numpy \
        python3-setuptools \
        python3-pip \
        rna-star \
        software-properties-common \
        tar \
        tzdata \
        unzip \
        vim \
        wget \
        zip \
        zlib1g \
        zlib1g-dev \
        && \
    rm -rf /var/lib/apt/lists/*

#### Descargar e instalar FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod +x FastQC/fastqc && \
    ln -s /FastQC/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.11.9.zip

#### Descargar e instalar MultiQC
RUN pip3 install --upgrade numpy
RUN pip3 install multiqc 

#### Descarga e instala Kallisto
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz && \
    tar -zxvf kallisto_linux-v0.46.1.tar.gz && \
    chmod +x kallisto/kallisto && \
    ln -s /kallisto/kallisto /usr/local/bin/kallisto && \
    rm kallisto_linux-v0.46.1.tar.gz

#### Descargar e instalar BWA
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar -xjvf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    ln -s /bwa-0.7.17/bwa /usr/local/bin/bwa && \
    cd .. && \
    rm bwa-0.7.17.tar.bz2

#### Descargar e instalar STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz && \
    tar -zxvf 2.7.9a.tar.gz && \
    cd STAR-2.7.9a/source && \
    make STAR && \
    ln -s /STAR-2.7.9a/source/STAR /usr/local/bin/STAR && \
    cd ../.. && \
    rm 2.7.9a.tar.gz

#### Descargar e Instalar GATK4
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip && \
    unzip gatk-4.2.6.1.zip && \
    chmod +x gatk-4.2.6.1/gatk && \
    ln -s /gatk-4.2.6.1/gatk /usr/local/bin/gatk && \
    rm gatk-4.2.6.1.zip

#### Descargar e instalar Picard    
RUN wget -q https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar && \
    mv picard.jar /usr/local/bin/ && \
    chmod +x /usr/local/bin/picard.jar

#### Descargar e instalar Trim Galore
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.7.tar.gz && \
    tar xvzf 0.6.7.tar.gz && \
    mv TrimGalore-0.6.7/trim_galore /usr/local/bin/ && \
    chmod +x /usr/local/bin/trim_galore && \
    rm -rf 0.6.7.tar.gz TrimGalore-0.6.7
    
#### Descargar e instalar Annovar 
#### IMPORTANTE: Debes de obtener tu propia liga para uso de investigación
#RUN wget link/to/annovar.latest.tar.gz && \
#    tar -zxvf annovar.latest.tar.gz && \
#    rm annovar.latest.tar.gz

#### Descargar e instalar Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 && \
    tar -xjvf samtools-1.12.tar.bz2 && \
    cd samtools-1.12 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm samtools-1.12.tar.bz2

#### Descargar e instalar Bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2 && \
    tar -xjvf bcftools-1.12.tar.bz2 && \
    cd bcftools-1.12 && \
    make && \
    make install && \
    cd .. && \
    rm bcftools-1.12.tar.bz2

### Instalar R version 4.2.3
RUN wget https://cdn.posit.co/r/ubuntu-2004/pkgs/r-4.2.3_1_amd64.deb && \
    apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -f -y ./r-4.2.3_1_amd64.deb && \
    ln -s /opt/R/4.2.3/bin/R /usr/bin/R && \
    ln -s /opt/R/4.2.3/bin/Rscript /usr/bin/Rscript && \
    ln -s /opt/R/4.2.3/lib/R /usr/lib/R && \
    rm r-4.2.3_1_amd64.deb && \
    rm -rf /var/lib/apt/lists/*

#### Instalar las librerías de R necesarias para los pipelines
RUN R -e "install.packages(c('data.table', 'dplyr', 'getopt', 'ggplot2', 'gplots', 'gsalib', 'optparse', 'backports', 'readr', 'ggrepel'), repos='http://cran.us.r-project.org')"
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
   libfftw3-dev \
   gcc && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
RUN R -e "install.packages('BiocManager',repos='http://cran.us.r-project.org')" -e "BiocManager::install(c('GenomicFeatures', 'GenomicAlignments', 'VariantAnnotation', 'S4Vectors', 'SummarizedExperiment', 'Biobase','tximportData', 'tximport', 'BUSpaRse', 'DESeq2', 'PCAtools', 'ComplexHeatmap', 'EnhancedVolcano', 'AnnotationDbi', 'topGO', 'GeneTonic', 'SummarizedExperiment', 'GenomicAlignments', 'org.Hs.eg.db', 'rhdf5', 'pcaExplorer'))"

RUN ln -s /usr/local/bin/picard.jar /usr/bin/picard.jar
#RUN ln -s /annovar/table_annovar.pl /usr/local/bin/table_annovar.pl
#RUN ln -s /annovar/table_annovar.pl /usr/bin/table_annovar.pl

### Configurar la variable de entorno para Java 8 ###
ENV JAVA_8_HOME "/usr/lib/jvm/java-8-openjdk-amd64"
ENV PATH "$JAVA_8_HOME/bin:$PATH"
ENV TZ "America/Mexico_City"

## Establecer el directorio de trabajo
WORKDIR /usr/local/lib/

# Ejecutar 
CMD ["/bin/bash"]
