#!/bin/sh
# Directorio de salida de nextflow
path=$1

## Ejecutar nextflow
nextflow run main.nf -resume -with-trace trace_vcs.txt -with-report report_vcs.html -with-timeline timeline_vcs.html

mkdir -p $path/run_files/modules

## Copiar los ejecutables, la configuración y los modulos de nextflow utilizados
cp ../../modules/VC-Somatic/vc-scommon/* $path/run_files/modules
cp ../../modules/VC-Somatic/vc-paired/* $path/run_files/modules
cp main.nf nextflow.config sample_info.tsv $path/run_files

## Mover los reportes de nextflow a una carpeta en el directrio de salida de nextflow
mv trace_vcs.txt timeline_vcs.html report_vcs.html $path/run_files

echo -e "Archivos y reportes del flujo de trabajo copiados a $path/run_files \n"
echo -e "Para explorar los resultados ejecuta: cd $path/out \n"
echo -e "Flujo de trabajo completado"
