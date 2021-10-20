#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.replicates = 1000
params.results = "results"

phylox = Channel.of("phylotorch", "phylotorch-bito", "phylojax")

process RUN_PHYSHER_MICRO {
  publishDir "$params.results/micro/physher.${size}"

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file)
  output:
  path("out.txt")
  """
  physher-benchmark ${params.replicates} ${seq_file} ${lsd_newick} > out.txt
  """
}

process RUN_PHYLOX_MICRO {
  publishDir "$params.results/micro/${phylox}.${size}.${rep}", mode: 'copy'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), val(phylox)
  output:
  path("out.txt")
  """
  python $baseDir/scripts/${phylox}/benchmark.py -i $seq_file \
                                                 -t $lsd_newick \
                                                 --replicates ${params.replicates} > out.txt
  """
}

workflow micro {
  take:
  data
  main:
  RUN_PHYSHER_MICRO(data)

//  RUN_PHYLOX_MICRO(data.combine(phylox))
}