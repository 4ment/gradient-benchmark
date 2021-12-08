#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.results = "results"
params.iterations = 5000

use_bito = Channel.of(true, false)


flu_H3N2 = "$baseDir/flu_H3N2"
physher_jc69_template = "$flu_H3N2/physher-JC69.template"
torchtree_jc69_template = "$flu_H3N2/phylotorch-JC69.template"


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
  publishDir "$params.results/macro/phylostan", mode: 'copy'

  input:
  tuple val(size), val(rep), path(tree_file), path(seq_file)
  path phylostan_stan
  path phylostan_pkl
  output:
  path 'phylostan', emit: phylostan_out
  path("out.txt")
  path("phylostan.${size}.${rep}.log")
  """
  { time \
  phylostan run -i ${seq_file} \
                -t ${tree_file} \
                -s ${phylostan_stan} \
                -o phylostan \
                -m JC69 \
                --heterochronous \
                --estimate_rate \
                --clock strict \
                 -c constant \
                --iter ${params.iterations}  \
                --eta 0.01 \
                --tol_rel_obj 0.00000001 \
                --elbo_samples 1 \
                --samples 1 > out.txt ; } 2> phylostan.${size}.${rep}.log & exit 0
  """
}

process PREPARE_PHYSHER {

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), path(lsd_dates)
  output:
  tuple val(size), val(rep), path("physher.json")
  """
  helper.py 1 \
                                    $seq_file \
                                    $lsd_newick \
                                    ${lsd_dates} \
                                    $physher_jc69_template physher.json \
                                    ${params.iterations}
  """
}

process RUN_PHYSHER {
  publishDir "$params.results/macro/physher", mode: 'copy'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), path(physher_json)
  output:
  path("out.txt")
  path("physher.${size}.${rep}.log")
  """
  { time physher $physher_json > out.txt ; } 2> physher.${size}.${rep}.log
  """
}

process PREPARE_TORCHTREE {
  label 'bito'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), path(lsd_dates), val(bito)
  output:
  tuple val(size), val(rep), path("torchtree.json"), val(bito)
  """
  helper.py 3 \
            $seq_file \
            $lsd_newick \
            ${lsd_dates} \
            $torchtree_jc69_template torchtree.json \
            ${params.iterations} \
            ${bito}
  """
}

process RUN_TORCHTREE {
  label 'bito'

  publishDir "$params.results/macro/torchtree", mode: 'copy'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), path(torchtree_json), val(bito)
  output:
  path("out.txt")
  path("torchtree.${bito}.${size}.${rep}.log")
  """
  { time \
  torchtree $torchtree_json > out.txt ; } 2> torchtree.${bito}.${size}.${rep}.log
  """
}

process RUN_PHYLOJAX {
  label 'bito'

  publishDir "$params.results/macro/phylojax", mode: 'copy'

  input:
  tuple val(size), val(rep), path(tree_file), path(seq_file)
  output:
  path("out.txt")
  path("phylojax.${size}.${rep}.log")
  """
  { time \
  phylojax -i ${seq_file} \
           -t ${tree_file} \
           --iter ${params.iterations} \
           --eta 0.01 \
           --elbo_samples 1 \
           --grad_samples 1 > out.txt ; } 2> phylojax.${size}.${rep}.log
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
}
