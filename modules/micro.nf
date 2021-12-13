#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.replicates = 1000
params.results = "results"

phylox = Channel.of("torchtree", "bitorch", "phylojax")

process RUN_PHYSHER_BENCHMARK {
  publishDir "$params.results/micro/physher", mode: 'copy'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file)
  output:
  path("physher.${size}.${rep}.csv")
  """
  physher-benchmark -i ${seq_file} \
                    -t ${lsd_newick} \
                    -r ${params.replicates} \
                    -s 0.001 \
                    -o out.csv
  awk 'NR==1{print "program,size,rep,precision,"\$0}; \
       NR>1{print "physher,$size,$rep,64,"\$0}' out.csv \
    > physher.${size}.${rep}.csv
  """
}

process RUN_PHYLOX_BENCHMARK {
  label 'bito'

  publishDir "$params.results/micro/${phylox}", mode: 'copy'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), val(phylox), val(precision)
  output:
  path("${phylox}.${size}.${rep}.${precision}.csv")
  script:
  if (precision == "32")
    extra = " -d float32"
  else
    extra = ""
  """
  ${phylox}-benchmark -i $seq_file \
                      -t $lsd_newick \
                      -r ${params.replicates} \
                      -s 0.001 \
                      -o out.csv \
                      ${extra}
  awk 'NR==1{print "program,size,rep,precision,"\$0}; \
       NR>1{print "$phylox,$size,$rep,${precision},"\$0}' out.csv \
      > ${phylox}.${size}.${rep}.${precision}.csv
  """
}

process COMBIME_CSV {
  publishDir "$params.results/micro/", mode: 'copy'

  input:
  path files
  output:
  path("micro.csv")

  """
  head -n1 ${files[0]} > micro
  tail -q -n+2 *.csv >> micro
  mv micro micro.csv
  """
}

workflow micro {
  take:
  data
  main:
  RUN_PHYSHER_BENCHMARK(data)

  RUN_PHYLOX_BENCHMARK(data.combine(phylox).combine(
          Channel.of("64")).mix(data.combine(Channel.of(['torchtree', "32"]))))

  ch_files = Channel.empty()
  ch_files = ch_files.mix(
          RUN_PHYSHER_BENCHMARK.out.collect(),
          RUN_PHYLOX_BENCHMARK.out.collect())
  COMBIME_CSV(ch_files.collect())
}
