FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

LABEL maintainer="Servicios INMEGEN"
LABEL version="1.0"

# Building takes several hours! Be sure to symlink /var/lib/docker to a larger disk.
# Build the image: docker build -t pipelinesinmegen/pipelines_inmegen:public -f Dockerfile .
# To run the container: docker run -i -t pipelinesinmegen/pipelines_inmegen:public /bin/bash

# Install linux tools
RUN apt-get update && \
    apt-get install -y \
        ant \
        automake \
        apt-utils \
        build-essential \
        bzip2 \
        ca-certificates \
        cutadapt \
        curl \
        gawk \
        gcc \
        git \
        g++ \
        hmmer \
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
        tabix \
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

#### Download and set FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    chmod a+x FastQC/fastqc && \
    ln -s /FastQC/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.12.1.zip

#### Install MultiQC
RUN pip3 install --upgrade numpy
RUN pip3 install multiqc 

#### Download and set Kallisto
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz && \
    tar -zxvf kallisto_linux-v0.46.1.tar.gz && \
    chmod a+x kallisto/kallisto && \
    ln -s /kallisto/kallisto /usr/local/bin/kallisto && \
    rm kallisto_linux-v0.46.1.tar.gz

#### Download and install BWA
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar -xjvf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    ln -s /bwa-0.7.17/bwa /usr/local/bin/bwa && \
    cd .. && \
    rm bwa-0.7.17.tar.bz2

#### Download and install STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz && \
    tar -zxvf 2.7.9a.tar.gz && \
    cd STAR-2.7.9a/source && \
    make STAR && \
    ln -s /STAR-2.7.9a/source/STAR /usr/local/bin/STAR && \
    cd ../.. && \
    rm 2.7.9a.tar.gz

#### Download and set GATK4
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip && \
    unzip gatk-4.2.6.1.zip && \
    chmod a+x gatk-4.2.6.1/gatk && \
    ln -s /gatk-4.2.6.1/gatk /usr/local/bin/gatk && \
    rm gatk-4.2.6.1.zip

#### Download and set Picard    
RUN wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar && \
    mkdir -p /picard/bin && \
    mv picard.jar /picard/ && \
    echo '#!/usr/bin/env bash' > /picard/bin/picard && \
    echo 'java -jar /picard/picard.jar ${@}' >> /picard/bin/picard && \
    chmod a+x /picard/bin/*

#### Download and set Trimmomatic    
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    mv Trimmomatic-0.39 trimmomatic && \
    cd /trimmomatic && \
    mv trimmomatic-0.39.jar trimmomatic.jar && \
    mkdir -p bin/ && \
    echo '#!/usr/bin/env bash' > /trimmomatic/bin/trimmomatic && \
    echo 'java -jar /trimmomatic/trimmomatic.jar ${@}' >> /trimmomatic/bin/trimmomatic && \
    chmod a+x /trimmomatic/bin/* && \
    cd ../ && \
    rm Trimmomatic-0.39.zip

#### Download and set Trim Galore
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.7.tar.gz && \
    tar xvzf 0.6.7.tar.gz && \
    mv TrimGalore-0.6.7/trim_galore /usr/local/bin/ && \
    chmod a+x /usr/local/bin/trim_galore && \
    rm -rf 0.6.7.tar.gz TrimGalore-0.6.7
    
#### Download and install Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 && \
    tar -xjvf samtools-1.12.tar.bz2 && \
    cd samtools-1.12 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm samtools-1.12.tar.bz2

#### Download and install Bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 && \
    tar -xjvf bcftools-1.19.tar.bz2 && \
    cd bcftools-1.19 && \
    make && \
    make install && \
    cd .. && \
    rm bcftools-1.19.tar.bz2

#### Download and install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 && \
    tar -xjvf htslib-1.19.1.tar.bz2 && \
    cd htslib-1.19.1 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm htslib-1.19.1.tar.bz2
 
#### Download and set mosdepth
RUN wget "https://github.com/brentp/mosdepth/releases/download/v0.3.6/mosdepth" && \
    cp ./mosdepth /usr/local/bin/ && \
    chmod a+x /usr/local/bin/mosdepth && \
    ln -s /usr/local/bin/mosdepth /usr/bin/mosdepth

#### Download and set FastQSreen
RUN wget https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.15.3.zip && \
    unzip v0.15.3.zip  && \
    chmod a+x FastQ-Screen-0.15.3/fastq_screen && \
    ln -s /FastQ-Screen-0.15.3/fastq_screen /usr/local/bin/fastq_screen && \
    rm v0.15.3.zip

#### Download and set Bowtie
RUN wget https://github.com/BenLangmead/bowtie/releases/download/v1.3.1/bowtie-1.3.1-linux-x86_64.zip  && \
    unzip bowtie-1.3.1-linux-x86_64.zip && \
    chmod a+x bowtie-1.3.1-linux-x86_64/bowtie && \
    ln -s /bowtie-1.3.1-linux-x86_64/bowtie /usr/local/bin/bowtie && \
    rm bowtie-1.3.1-linux-x86_64.zip

#### Download and set Bowtie2
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.3/bowtie2-2.5.3-linux-x86_64.zip  && \
    unzip bowtie2-2.5.3-linux-x86_64.zip && \
    chmod a+x bowtie2-2.5.3-linux-x86_64/bowtie2 && \
    ln -s /bowtie2-2.5.3-linux-x86_64/bowtie2 /usr/local/bin/bowtie2 && \
    rm bowtie2-2.5.3-linux-x86_64.zip

#### Download and set Fastp
RUN wget http://opengene.org/fastp/fastp && \
    cp ./fastp /usr/local/bin/ && \
    chmod a+x /usr/local/bin/fastp && \
    ln -s /usr/local/bin/fastp /usr/bin/fastp

#### Download and set Mothur
RUN wget https://github.com/mothur/mothur/releases/download/v1.48.0/Mothur.Ubuntu_20.zip && \
    unzip Mothur.Ubuntu_20.zip && \
    chmod a+x ./mothur/mothur && \
    ln -s /mothur/mothur /usr/local/bin/mothur && \
    ln -s /mothur/mothur /usr/bin/mothur && \
    rm Mothur.Ubuntu_20.zip

#### Symbolic link
RUN ln -s /usr/bin/python3 /usr/bin/python

### Download and install R 
RUN wget https://cdn.posit.co/r/ubuntu-2004/pkgs/r-4.3.2_1_amd64.deb && \
    apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -f -y ./r-4.3.2_1_amd64.deb && \
    ln -s /opt/R/4.3.2/bin/R /usr/bin/R && \
    ln -s /opt/R/4.3.2/bin/Rscript /usr/bin/Rscript && \
    ln -s /opt/R/4.3.2/lib/R /usr/lib/R && \
    rm r-4.3.2_1_amd64.deb && \
    rm -rf /var/lib/apt/lists/*

#### Install R packages
RUN R -e "install.packages(c('devtools','data.table','tidyr','dplyr','getopt','gplots','gsalib','optparse','backports','readr','ggrepel','mutSignatures','enrichR','vroom'), repos='http://cran.us.r-project.org')"
RUN apt-get update && apt-get install -y --no-install-recommends \
   libfftw3-dev \
   gcc && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
RUN R -e "install.packages('BiocManager',repos='http://cran.us.r-project.org')" -e "BiocManager::install(c('GenomicFeatures', 'GenomicAlignments', 'VariantAnnotation', 'S4Vectors', 'SummarizedExperiment', 'Biobase','tximportData', 'tximport', 'BUSpaRse', 'DESeq2', 'edgeR', 'dada2', 'PCAtools', 'ComplexHeatmap', 'maftools' , 'pheatmap', 'NMF', 'EnhancedVolcano', 'AnnotationDbi', 'topGO', 'GeneTonic', 'SummarizedExperiment', 'GenomicAlignments', 'BiocGenerics', 'org.Hs.eg.db', 'rhdf5', 'pcaExplorer','BSgenome.Hsapiens.UCSC.hg38','IRanges','GenomicRanges','BSgenome','GenomeInfoDb'))"
RUN R -e "cat(\"local(options(repos = c(CRAN = 'https://cloud.r-project.org')))\n\", file = R.home('etc/Rprofile.site'), append = TRUE)"
RUN R -e "install.packages('remotes')"
RUN R -e "remotes::install_version('yaml', version = '2.3.7', upgrade = FALSE)"
RUN R -e "remotes::install_version('vroom', version = '1.6.5', upgrade = FALSE)"
RUN R -e "remotes::install_version('plotly', version = '4.10.3', upgrade = FALSE)"
RUN R -e "remotes::install_version('ggplot2', version = '3.4.4', upgrade = FALSE)"
Run R -e "remotes::install_version('scales', version = '1.2.1', upgrade = FALSE)"
RUN R -e "BiocManager::install(c('cn.mops'))"
RUN R -e "install.packages(c('checkmate', 'foreign', 'htmlTable', 'viridis', 'Formula', 'Hmisc'), repos='http://cran.r-project.org')"

### Environment variables
ENV JAVA_8_HOME "/usr/lib/jvm/java-8-openjdk-amd64"
ENV PATH "$JAVA_8_HOME/bin:$PATH"
ENV PATH /picard/bin:/trimmomatic/bin:${PATH}

ENV TZ "America/Mexico_City"

# Set run command 
CMD ["/bin/bash"]
