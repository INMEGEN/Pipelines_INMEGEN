process starcmatrix {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  containerOptions "-v ${params.refdir}:/ref"
  publishDir params.out + "/StarMCounts" + params.DEAname , mode: 'copy'

  input:
  val(samples_s)
  path(star_dir)

  output:
  path("star_mcounts.tsv"), emit: mcounts


  script:
  """
   > "star_mcounts.tsv"

  output_file="star_mcounts_tmp.tsv"

  first_file=\$(ls ${star_dir}/*.tab | head -n 1) 
  awk -F'\t' '{print \$1, \$2}' "\$first_file" > "\$output_file"

  for file in ${star_dir}/*.tab; do
    if [ "\$file" != "\$first_file" ]; then
    paste -d '\t' "\$output_file" <(awk -F'\t' '{print \$2}' "\$file") > temp && mv temp "\$output_file"
    fi
  done

  awk '{\$1=\$1}1' OFS="\t" star_mcounts_tmp.tsv > star_mcounts.tsv

  rm star_mcounts_tmp.tsv
  """
}
