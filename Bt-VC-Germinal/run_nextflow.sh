#!/bin/sh
# Directorio de salida de nextflow
path=$1

## Ejecutar nextflow
nextflow run main.nf -resume -with-trace trace_VCG.txt -with-report report_VCG.html -with-timeline timeline_VCG.html

mkdir -p $path/run_files/modules

## Copiar los ejecutables, la configuración y los modulos de nextflow utilizados
cp ../modules/VC-Germinal/* $path/run_files/modules
cp ../modules/qualitycontrol/* $path/run_files/modules
cp ../modules/bootstrapping/* $path/run_files/modules
cp ../modules/metricts/* $path/run_files/modules
cp ../modules/common/* $path/run_files/modules
cp ../modules/annotation/* $path/run_files/modules
cp main.nf nextflow.config sample_info.tsv $path/run_files

## Mover los reportes de nextflow a una carpeta en el directrio de salida de nextflow
mv trace_VCG.txt timeline_VCG.html report_VCG.html $path/run_files

echo -e "Archivos y reportes del flujo de trabajo copiados a $path/run_files \n"
echo -e "Para explorar los resultados ejecuta: cd $path/out \n"
echo -e "Flujo de trabajo completado"
