process salmon_version {
  cache 'lenient'
  container 'combinelab/salmon:latest'
  publishDir params.out +"/salmon", mode: 'copy'

  input:
  file (int1)

  output:
  path("salmon_version.txt"), emit: sversion

  script:
  """
  echo "\$(salmon --version)" > salmon_version.txt 
  """
}
