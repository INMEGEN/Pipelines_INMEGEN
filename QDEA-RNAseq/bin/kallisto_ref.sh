#!/bin/bash

### Los siguientes comandos generan el índice de kallisto, descargan el archivo gtf y generan la referencia de star.

mkdir kallisto_idx/

cd kallisto_idx/

wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

cd ..

docker run --cpus 6 --user="$(id -u):$(id -g)" -v "$PWD/kallisto_idx":/data pipelinesinmegen/pipelines_inmegen:public kallisto index -i /data/Homo_sapiens.GRCh38.idx /data/Homo_sapiens.GRCh38.cdna.all.fa.gz

mkdir star_ref/

cd star_ref/

https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz

cd ..

docker run --cpus 10 -v "$PWD/star_ref":/data pipelinesinmegen/pipelines_inmegen:public STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data --genomeFastaFiles /data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --sjdbGTFfile /data/Homo_sapiens.GRCh38.111.gtf.gz

## Recuerda modificar la propiedad de la referencia de star
