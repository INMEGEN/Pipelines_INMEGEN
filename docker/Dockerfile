FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

LABEL maintainer="Servicios INMEGEN"
LABEL version="1.0"

# Building takes several hours! Be sure to symlink /var/lib/docker to a larger disk.
# Build the image: docker build -t pipelines_inmegen:public -f Dockerfile .
# To run the container: docker run -i -t pipelines_inmegen:public /bin/bash

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
RUN wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar && \
    mv picard.jar /usr/local/bin/ && \
    chmod +x /usr/local/bin/picard.jar

#### Descargar e instalar Trimmomatic    
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    cd Trimmomatic-0.39 && \
    mv trimmomatic-0.39.jar /usr/local/bin/ && \
    chmod +x /usr/local/bin/trimmomatic-0.39.jar && \
    cd ../ && \
    rm Trimmomatic-0.39.zip

#### Descargar e instalar Trim Galore
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.7.tar.gz && \
    tar xvzf 0.6.7.tar.gz && \
    mv TrimGalore-0.6.7/trim_galore /usr/local/bin/ && \
    chmod +x /usr/local/bin/trim_galore && \
    rm -rf 0.6.7.tar.gz TrimGalore-0.6.7
    
#### Descargar e instalar Annovar
##RUN wget -URL to annovar.latest.tar.gz- && \
##    tar -zxvf annovar.latest.tar.gz && \
##    rm annovar.latest.tar.gz

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
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 && \
    tar -xjvf bcftools-1.19.tar.bz2 && \
    cd bcftools-1.19 && \
    make && \
    make install && \
    cd .. && \
    rm bcftools-1.19.tar.bz2

#### Descargar e instalar htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 && \
    tar -xjvf htslib-1.19.1.tar.bz2 && \
    cd htslib-1.19.1 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm htslib-1.19.1.tar.bz2
 
#### Descargar e instalar mosdepth
RUN wget "https://github.com/brentp/mosdepth/releases/download/v0.3.6/mosdepth" && \
    cp ./mosdepth /usr/local/bin/ && \
    ln -s /usr/local/bin/mosdepth /usr/bin/mosdepth && \
    chmod +x /usr/local/bin/mosdepth

### Instalar R 
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

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
   cutadapt && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

#### Descargar e instalar FastQSreen
RUN wget https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.15.3.zip && \
    unzip v0.15.3.zip  && \
    chmod +x FastQ-Screen-0.15.3/fastq_screen && \
    ln -s /FastQ-Screen-0.15.3/fastq_screen /usr/local/bin/fastq_screen && \
    rm v0.15.3.zip

#### Descargar e instalar Bowtie
RUN wget https://github.com/BenLangmead/bowtie/releases/download/v1.3.1/bowtie-1.3.1-linux-x86_64.zip  && \
    unzip bowtie-1.3.1-linux-x86_64.zip && \
    chmod +x bowtie-1.3.1-linux-x86_64/bowtie && \
    ln -s /bowtie-1.3.1-linux-x86_64/bowtie /usr/local/bin/bowtie && \
    rm bowtie-1.3.1-linux-x86_64.zip

#### Descargar e instalar Bowtie2
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.3/bowtie2-2.5.3-linux-x86_64.zip  && \ 
    unzip bowtie2-2.5.3-linux-x86_64.zip && \
    chmod +x bowtie2-2.5.3-linux-x86_64/bowtie2 && \
    ln -s /bowtie2-2.5.3-linux-x86_64/bowtie2 /usr/local/bin/bowtie2 && \ 
    rm bowtie2-2.5.3-linux-x86_64.zip

#RUN chmod 774 /usr/local/bin/picard.jar
RUN ln -s /usr/local/bin/picard.jar /usr/bin/picard.jar
RUN ln -s /usr/local/bin/trimmomatic-0.39.jar /usr/bin/trimmomatic-0.39.jar
RUN ln -s /annovar/table_annovar.pl /usr/local/bin/table_annovar.pl
RUN ln -s /annovar/table_annovar.pl /usr/bin/table_annovar.pl
RUN ln -s /usr/bin/python3 /usr/bin/python

### Configurar la variable de entorno para Java 8 ###
ENV JAVA_8_HOME "/usr/lib/jvm/java-8-openjdk-amd64"
ENV PATH "$JAVA_8_HOME/bin:$PATH"
ENV TZ "America/Mexico_City"

## Establecer el directorio de trabajo
WORKDIR /usr/local/lib/

#ADD /usr/local/bin/picard.jar /usr/local/lib/

# Ejecutar 
CMD ["/bin/bash"]
