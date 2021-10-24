#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.replicates = 1000
params.results = "results"

phylox = Channel.of("phylotorch", "bitorch", "phylojax")

process RUN_PHYSHER_BENCHMARK {
  publishDir "$params.results/micro/physher", mode: 'copy'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file)
  output:
  path("physher.${size}.${rep}.csv")
  """
  physher-benchmark ${params.replicates} ${seq_file} ${lsd_newick} \
    | physher-parser.py - \
    | awk 'NR==1{print "program,size,rep,"\$0};NR>1{print "physher,$size,$rep,"\$0}' \
    > physher.${size}.${rep}.csv
  """
}

process RUN_PHYLOX_BENCHMARK {
  publishDir "$params.results/micro/${phylox}", mode: 'copy'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), val(phylox)
  output:
  path("${phylox}.${size}.${rep}.csv")
  """
  source activate bito
  ${phylox}-benchmark -i $seq_file \
                      -t $lsd_newick \
                      --replicates ${params.replicates} \
                      -o out.csv
  awk 'NR==1{print "program,size,rep,"\$0};NR>1{print "$phylox,$size,$rep,"\$0}' out.csv \
    > ${phylox}.${size}.${rep}.csv
  """
}

process COMBIME_CSV {
  publishDir "$params.results/micro/", mode: 'copy'

  input:
  path files
  output:
  path("micro.csv")

  """
  head -n1 ${files[0]} > micro.csv
  tail -q -n+2 *.csv >> micro.csv
  """
}

workflow micro {
  take:
  data
  main:
  RUN_PHYSHER_BENCHMARK(data)

  RUN_PHYLOX_BENCHMARK(data.combine(phylox))

  ch_files = Channel.empty()
  ch_files = ch_files.mix(
          RUN_PHYSHER_BENCHMARK.out.collect(),
          RUN_PHYLOX_BENCHMARK.out.collect())
  COMBIME_CSV(ch_files.collect())
}