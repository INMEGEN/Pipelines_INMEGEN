#!/bin/bash

### Usa los siguientes comandos para gennerar el índice de kallisto y descargar el archivo gtf

mkdir kallisto_ref

cd kallisto_ref/

wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

kallisto index -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx Homo_sapiens.GRCh38.cdna.all.fa.gz

wget https://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
