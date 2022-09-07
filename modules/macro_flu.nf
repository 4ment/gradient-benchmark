#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.results = "results"
params.iterations = 5000

include { COMBIME_TIME_LOG } from './utils.nf'

use_bito = Channel.of(true, false)

flu_H3N2 = "$baseDir/flu_H3N2"
physher_jc69_template = "$flu_H3N2/physher-JC69.template"


process COMPILE_PHYLOSTAN {
  label 'fast'

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
                  --clockpr exponential \
                  -c constant \
                  --compile
  """
}

process RUN_PHYLOSTAN {
  errorStrategy 'ignore'

  publishDir "$params.results/macro/phylostan", mode: 'copy'

  input:
    tuple val(size), val(rep), path(tree_file), val(rate), path(seq_file)
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
                --heights_init \
                --rate ${rate} \
                --samples 1 > phylostan.${size}.${rep}.txt ; } 2> phylostan.${size}.${rep}.log
  """
}

process PREPARE_PHYSHER {
  label 'ultrafast'

  publishDir "$params.results/macro/physher", mode: 'copy'

  input:
    tuple val(size), val(rep), path(lsd_newick), val(rate), path(seq_file), path(lsd_dates)
  output:
    tuple val(size), val(rep), path("physher.${size}.${rep}.json")
  """
  helper.py physher --input $seq_file \
                    --tree $lsd_newick \
                    --dates ${lsd_dates} \
                    --template $physher_jc69_template \
                    --output physher.${size}.${rep}.json \
                    --iterations ${params.iterations} \
                    --rate ${rate} \
                    --lr 0.0001 \
                    --tol 0.0000001
  """
}

process RUN_PHYSHER {
  label 'fast'

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
  label 'ultrafast'
  label 'bito'

  publishDir "$params.results/macro/torchtree", mode: 'copy'

  input:
    tuple val(size), val(rep), path(lsd_newick), val(rate), path(seq_file), path(lsd_dates), val(bito)
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
                --rate_init ${rate} \
                --clockpr exponential \
                --lr 0.0001 \
                --elbo_samples 1 \
                --tol_rel_obj 0 \
                --iter ${params.iterations} \
                --samples 0 \
                ${bito_arg} > torchtree.${bito}.${size}.${rep}.json
  """
}

process RUN_TORCHTREE {
  label 'normal'
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
    tuple val(size), val(rep), path(tree_file), val(rate), path(seq_file)
  output:
    tuple path("phylojax.${size}.${rep}.txt"), path("phylojax.${size}.${rep}.log")

  when:
    size <= 750
  """
  { time \
  phylojax -i ${seq_file} \
           -t ${tree_file} \
           --iter ${params.iterations} \
           --eta 0.01 \
           --elbo_samples 1 \
           --rate_init ${rate} \
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

workflow macro_flu {
  take:
  data
  main:
  COMPILE_PHYLOSTAN("H3N2_HA_2011_2013", "JC69")

  // this channel does not contain the substitution rate
  data_run = data.map { tuple(it[0], it[1], it[2], it[4]) }

  RUN_PHYLOSTAN(
          data.map { it.take(5) },
          COMPILE_PHYLOSTAN.out.stan,
          COMPILE_PHYLOSTAN.out.pkl)

  PREPARE_PHYSHER(data)

  RUN_PHYSHER(data_run.join(PREPARE_PHYSHER.out, by: [0, 1]))

  RUN_PHYLOJAX(data.map { it.take(5) })

  PREPARE_TORCHTREE(data.combine(use_bito))

  RUN_TORCHTREE(data_run.join(PREPARE_TORCHTREE.out, by: [0, 1]))

  RUN_TREEFLOW(data_run)

  ch_files = Channel.empty()
  ch_files = ch_files.mix(
        RUN_PHYSHER.out.collect(),
        RUN_PHYLOJAX.out.collect(),
        RUN_PHYLOSTAN.out[0].collect(),
        RUN_TORCHTREE.out.collect(),
        RUN_TREEFLOW.out.collect())

  COMBIME_TIME_LOG(ch_files.collect(), "macro")
}
