#!/bin/bash

### Los siguientes comandos generan el índice de kallisto y descargan el archivo gtf.

mkdir hg38/

cd hg38/

wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

wget https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz

cd ..

docker run --cpus 6 -v "$PWD/hg38":/data pipelinesinmegen/pipelines_inmegen:public kallisto index -i /data/Homo_sapiens.GRCh38.cdna.all.idx /data/Homo_sapiens.GRCh38.cdna.all.fa.gz

