#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.results = "results"
params.iterations = 5000

use_bito = Channel.of(true, false)

clock_rate = 0.002 // could be parsed from lsd ouput file

flu_H3N2 = "$baseDir/flu_H3N2"
physher_jc69_template = "$flu_H3N2/physher-JC69.template"


process COMPILE_PHYLOSTAN {
  input:
  val(name)
  val(model)
  output:
  path "${name}.stan", emit: stan
  path "${name}.pkl", emit: pkl
  """
  phylostan build -s ${name}.stan \
                  -m ${model} \
                  --heterochronous \
                  --estimate_rate \
                  --clock strict \
                  -c constant \
                  --compile
  """
}

process RUN_PHYLOSTAN {
  errorStrategy 'ignore'

  publishDir "$params.results/macro/phylostan", mode: 'copy'

  input:
  tuple val(size), val(rep), path(tree_file), path(seq_file)
  path phylostan_stan
  path phylostan_pkl
  output:
  tuple path("phylostan.${size}.${rep}.txt"), path("phylostan.${size}.${rep}.log")
  path "phylostan.${size}.${rep}", emit: phylostan_out
  path("phylostan.${size}.${rep}.diag")
  """
  { time \
  phylostan run -i ${seq_file} \
                -t ${tree_file} \
                -s ${phylostan_stan} \
                -o phylostan.${size}.${rep} \
                -m JC69 \
                --heterochronous \
                --estimate_rate \
                --clock strict \
                --clockpr exponential \
                -c constant \
                --iter ${params.iterations}  \
                --eta 0.0001 \
                --tol_rel_obj 0.00000001 \
                --elbo_samples 1 \
                --samples 1 > phylostan.${size}.${rep}.txt ; } 2> phylostan.${size}.${rep}.log
  """
}

process PREPARE_PHYSHER {
  publishDir "$params.results/macro/physher", mode: 'copy'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), path(lsd_dates)
  output:
  tuple val(size), val(rep), path("physher.${size}.${rep}.json")
  """
  helper.py 1 \
                                    $seq_file \
                                    $lsd_newick \
                                    ${lsd_dates} \
                                    $physher_jc69_template physher.${size}.${rep}.json \
                                    ${params.iterations} \
                                    ${clock_rate}
  """
}

process RUN_PHYSHER {
  publishDir "$params.results/macro/physher", mode: 'copy'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), path(physher_json)
  output:
  tuple path("physher.${size}.${rep}.txt"), path("physher.${size}.${rep}.log")
  """
  { time physher $physher_json > physher.${size}.${rep}.txt ; } 2> physher.${size}.${rep}.log
  """
}

process PREPARE_TORCHTREE {
  label 'bito'

  publishDir "$params.results/macro/torchtree", mode: 'copy'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), path(lsd_dates), val(bito)
  output:
  tuple val(size), val(rep), path("torchtree.${bito}.${size}.${rep}.json"), val(bito)
  script:
  if (bito)
    bito_arg = " --engine bitorch"
  else
    bito_arg = ""
  """
  torchtree-cli advi \
                -i $seq_file \
                -t $lsd_newick \
                --clock strict \
                --coalescent constant \
                --heights_init tree \
                --rate_init ${clock_rate} \
                --clockpr exponential \
                --eta 0.0001 \
                --elbo_samples 1 \
                --tol_rel_obj 0 \
                --iter ${params.iterations} \
                --samples 0 \
                ${bito_arg} > torchtree.${bito}.${size}.${rep}.json
  """
}

process RUN_TORCHTREE {
  label 'bito'

  errorStrategy 'ignore'

  publishDir "$params.results/macro/torchtree", mode: 'copy'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), path(torchtree_json), val(bito)
  output:
  tuple path("torchtree.${bito}.${size}.${rep}.txt"), path("torchtree.${bito}.${size}.${rep}.log")
  """
  { time \
  torchtree $torchtree_json > torchtree.${bito}.${size}.${rep}.txt ; } 2> torchtree.${bito}.${size}.${rep}.log
  """
}

process RUN_PHYLOJAX {
  label 'bito'

  errorStrategy 'ignore'

  publishDir "$params.results/macro/phylojax", mode: 'copy'

  input:
  tuple val(size), val(rep), path(tree_file), path(seq_file)
  output:
  tuple path("phylojax.${size}.${rep}.txt"), path("phylojax.${size}.${rep}.log")
  """
  { time \
  phylojax -i ${seq_file} \
           -t ${tree_file} \
           --iter ${params.iterations} \
           --eta 0.01 \
           --elbo_samples 1 \
           --rate_init ${clock_rate} \
           --heights_init tree \
           --grad_samples 1 > phylojax.${size}.${rep}.txt ; } 2> phylojax.${size}.${rep}.log
  """
}

process RUN_TREEFLOW {
  label 'bito'

  publishDir "$params.results/macro/treeflow", mode: 'copy'

  input:
  tuple val(size), val(rep), path(tree_file), path(seq_file)
  output:
  tuple path("treeflow.${size}.${rep}.txt"), path("treeflow.${size}.${rep}.log")
  """
  { time \
  treeflow_vi -i ${seq_file} \
              -t ${tree_file} \
              -n ${params.iterations} > treeflow.${size}.${rep}.txt ; } 2> treeflow.${size}.${rep}.log
  """
}

process COMBIME_TIME_LOG {
  publishDir "$params.results/macro/", mode: 'copy'

  input:
  path files
  output:
  path("macro.csv")

  """
  time-parser.py macro.csv ${files}
  """
}

workflow macro_flu {
  take:
  data
  main:
  COMPILE_PHYLOSTAN("H3N2_HA_2011_2013", "JC69")

  RUN_PHYLOSTAN(
          data.map { it.take(4) },
          COMPILE_PHYLOSTAN.out.stan,
          COMPILE_PHYLOSTAN.out.pkl)

  PREPARE_PHYSHER(data)

  RUN_PHYSHER(data.map { it.take(4) }.join(PREPARE_PHYSHER.out, by: [0, 1]))

  RUN_PHYLOJAX(data.map { it.take(4) })

  PREPARE_TORCHTREE(data.combine(use_bito))

  RUN_TORCHTREE(data.map { it.take(4) }.join(PREPARE_TORCHTREE.out, by: [0, 1]))

  RUN_TREEFLOW(data.map { it.take(4) })

  ch_files = Channel.empty()
  ch_files = ch_files.mix(
        RUN_PHYSHER.out.collect(),
        RUN_PHYLOJAX.out.collect(),
        RUN_PHYLOSTAN.out[0].collect(),
        RUN_TORCHTREE.out.collect(),
        RUN_TREEFLOW.out.collect())

  COMBIME_TIME_LOG(ch_files.collect())
}
