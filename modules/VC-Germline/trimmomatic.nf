process trimmomatic {
  cache 'lenient'
  publishDir params.out, mode: 'symlink'

  input:
  tuple val(sample), val(sample_id), val(PU), val(PL), val(LB), path(R1), path(R2)
  file(adapters)

  output:
  tuple val(sample), val(sample_id), val(PU), val(PL), val(LB), path("trimming_files/${sample_id}_R1.trimmed.fastq.gz"), path("trimming_files/${sample_id}_R2.trimmed.fastq.gz"), emit: trim_fq
  path("trimming_files/*un.trimmed.fastq.gz")

  script:
  """
  mkdir -p trimming_files
  cp ${R1} ${R2} trimming_files/

  docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v \$PWD/trimming_files:/data pipelinesinmegen/pipelines_inmegen:public \
  java -jar /usr/bin/trimmomatic-0.39.jar PE -threads ${params.ncrs} /data/${R1} /data/${R2} \
  /data/${sample_id}_R1.trimmed.fastq.gz /data/${sample_id}_R1un.trimmed.fastq.gz \
  /data/${sample_id}_R2.trimmed.fastq.gz /data/${sample_id}_R2un.trimmed.fastq.gz \
  ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30;

  cd trimming_files/
  rm ${R1} ${R2}
  """
}
