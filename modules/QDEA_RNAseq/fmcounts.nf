process fmcounts {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  publishDir params.out + "/fmcounts", mode: 'copy'

  input:
  val(samples)
  path(quants)

  output:
  path("fmcounts.tsv")      , emit: mcounts
  path("fmcounts_tpm.tsv")  , emit: tpm_mcounts

  script:
  """
  > "fmcounts_tmp.tsv"

  > "fmcounts_TPM_tmp.tsv"

  output_file="fmcounts_tmp.tsv"

  output_file_TPM="fmcounts_TPM_tmp.tsv"

  first_file=\$(ls ${quants}/*.tsv | head -n 1) 
  awk -F'\t' '{print \$1, \$3}' "\$first_file" > "\$output_file"

  for file in ${quants}/*.tsv; do
    if [ "\$file" != "\$first_file" ]; then
    paste -d '\t' "\$output_file" <(awk -F'\t' '{print \$3}' "\$file") >> temp && mv temp "\$output_file"
    fi
  done

  awk '{\$1=\$1}1' OFS='\t' fmcounts_tmp.tsv > fmcounts.tsv

  awk ' NR == 1 {sample_name = \$3; next} 
        NR == 2 {print "GeneID\\t" sample_name} 
        NR > 1 {gene_id = \$1; length_kb = \$2 / 1000; if (length_kb > 0) {rpk = \$3 / length_kb } else {rpk = 0}; sum_rpk += rpk; rpk_values[gene_id] = rpk} END {factor = (sum_rpk > 0) ? sum_rpk / 1000000 : 1; for (gene in rpk_values) { tpm = rpk_values[gene] / factor; print gene "\\t" tpm } 
     }' "\$first_file" > "\$output_file_TPM"  

  for file2 in ${quants}/*_counts.tsv; do
      if [ "\$file2" != "\$first_file" ]; then
      awk ' NR == 1 {sample_name = \$3; next}
            NR == 2 {print sample_name}
	    BEGIN {sum_rpk = 0;}
            NR > 1 { gene_id = \$1; length_kb = \$2 / 1000; if (length_kb > 0) {rpk = \$3 / length_kb } else {rpk = 0};  sum_rpk += rpk; rpk_values[gene_id] = rpk} END {factor = (sum_rpk > 0) ? sum_rpk / 1000000 : 1; for (gene in rpk_values) {tpm = rpk_values[gene] / factor; print tpm } 
         }' "\$file2" > temp_tpm.tsv
      paste -d '\t' "\$output_file_TPM" <(awk -F'\t' '{print \$1}' temp_tpm.tsv) > temp_tpm && mv temp_tpm "\$output_file_TPM"
      fi
  done 

  awk '{\$1=\$1}1' OFS='\t' fmcounts_TPM_tmp.tsv > fmcounts_tpm.tsv

  rm fmcounts_tmp.tsv fmcounts_TPM_tmp.tsv
  """
}
