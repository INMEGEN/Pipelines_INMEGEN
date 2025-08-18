#!/bin/bash

# Verificar que se pasó un argumento
if [[ -z "$1" ]]; then
    echo "Uso: $0 /ruta/al/directorio"
    exit 1
fi

DIR="$1"
OUTPUT="sample_info.tsv"

# Escribir encabezado
echo -e "Sample\tPath\tIndex" > "$OUTPUT"

# Buscar todos los .vcf.gz con su .tbi correspondiente
find "$DIR" -type f -name "*.vcf.gz" | while read -r vcf; do
    index="${vcf}.tbi"
    if [[ -f "$index" ]]; then
        sample=$(basename "$vcf" .vcf.gz)
        echo -e "${sample}\t${vcf}\t${index}" >> "$OUTPUT"
    fi
done
