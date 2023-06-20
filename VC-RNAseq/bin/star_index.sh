#!/usr/bin/env bash

mkdir -p ref/STAR_index
cd ref/STAR_index

#### Los archivos fasta y gtf se obtuvieron de  https://www.gencodegenes.org/human/
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz

gunzip gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

#### Crear diccionario y fasta index
samtools faidx GRCh38.primary_assembly.genome.fa
java -jar /path/to/picard.jar CreateSequenceDictionary -R GRCh38.primary_assembly.genome.fa -O GRCh38.primary_assembly.genome.dict

#### STAR index
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v43.chr_patch_hapl_scaff.annotation.gtf
