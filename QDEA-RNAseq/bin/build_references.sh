#!/bin/sh

mkdir -p references/

##### Download resources (dna, cdna and gtf file)
cd references

wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz & gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz & gunzip Homo_sapiens.GRCh38.cdna.all.fa.g

wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz & gunzip Homo_sapiens.GRCh38.113.gtf.gz

##### Reference STAR
mkdir -p STAR_index

cp Homo_sapiens.GRCh38.dna.primary_assembly.fa STAR_index/

cp Homo_sapiens.GRCh38.113.gtf STAR_index/

docker run --cpus 10  -v "$PWD/STAR_index":/data pipelinesinmegen/pipelines_inmegen:public STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data --genomeFastaFiles /data/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /data/Homo_sapiens.GRCh38.113.gtf --sjdbOverhang 99

#### Salmon index
mkdir -p Salmon_index

mv Homo_sapiens.GRCh38.dna.primary_assembly.fa Salmon_index/

mv Homo_sapiens.GRCh38.cdna.all.fa Salmon_index/

mv Homo_sapiens.GRCh38.113.gtf Salmon_index/

cd Salmon_index/

cat Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.dna.primary_assembly.fa > gentrome.fa

grep "^>" salmon_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

cd ..

docker run --cpus 10 -v "$PWD/Salmon_index":/data combinelab/salmon:latest salmon index -p 10 -t /data/gentrome.fa -d /data/decoys.txt -i /data/GRCh38_salmon_index -k 31
