#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.results = "results"
params.iterations = 5000

use_bito = Channel.of(true, false)


flu_H3N2 = "$baseDir/flu_H3N2"
physher_jc69_template = "$flu_H3N2/physher-JC69.template"
phylotorch_jc69_template = "$flu_H3N2/phylotorch-JC69.template"

params.subtrees_alignment = "$baseDir/treetime_validation/resources/flu_H3N2/H3N2_HA_2011_2013.fasta"


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
  publishDir "$params.results/phylostan.${size}", mode: 'copy'

  input:
  tuple val(size), val(rep), path(tree_file), path(seq_file)
  path phylostan_stan
  path phylostan_pkl
  output:
  path 'phylostan', emit: phylostan_out
  path("out.txt")
  path("time.log")
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
                --samples 1 > out.txt ; } 2> time.log
  """
}

process PREPARE_PHYSHER {

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), path(lsd_dates)
  output:
  tuple val(size), val(rep), path("physher.json")
  """
  python $baseDir/scripts/helper.py 1 \
                                    $seq_file \
                                    $lsd_newick \
                                    ${lsd_dates} \
                                    $physher_jc69_template physher.json \
                                    ${params.iterations}
  """
}

process RUN_PHYSHER {
  publishDir "$params.results/physher.${size}", mode: 'copy'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), path(physher_json)
  output:
  path("out.txt")
  path("time.log")
  """
  { time physher $physher_json > out.txt ; } 2> time.log
  """
}

process PREPARE_PHYLOTORCH {

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), path(lsd_dates), val(bito)
  output:
  tuple val(size), val(rep), path("phylotorch.json"), val(bito)
  """
  python $baseDir/scripts/helper.py 3 \
                                    $seq_file \
                                    $lsd_newick \
                                    ${lsd_dates} \
                                    $phylotorch_jc69_template phylotorch.json \
                                    ${params.iterations} \
                                    ${bito}
  """
}

process RUN_PHYLOTORCH {
  publishDir "$params.results/phylotorch.${bito}.${size}", mode: 'copy'

  input:
  tuple val(size), val(rep), path(lsd_newick), path(seq_file), path(phylotorch_json), val(bito)
  output:
  path("out.txt")
  path("time.log")
  """
  { time \
  phylotorch $phylotorch_json > out.txt ; } 2> time.log
  """
}

process RUN_PHYLOJAX {
  publishDir "$params.results/phylojax.${size}.${rep}", mode: 'copy'

  input:
  tuple val(size), val(rep), path(tree_file), path(seq_file)
  output:
  path("out.txt")
  path("time.log")
  """
  { time \
  phylojax -i ${seq_file} \
           -t ${tree_file} \
           --iter ${params.iterations} \
           --eta 0.01 \
           --elbo_samples 1 \
           --grad_samples 1 > out.txt ; } 2> time.log
  """
}

process RUN_LSD {
  input:
  tuple val(size), val(rep), path(subtree_dates_file), path(fasta), path(new_dates_file)
  output:
  tuple val(size), val(rep), path("lsd.out.date.nexus")
  path "lsd.out.nwk", emit: lsd_tree_newick // branch=subst
  path "lsd.out"

  """
  lsd2 -i ${subtree_dates_file} \
       -d ${new_dates_file} \
       -o lsd.out \
       -s 1701
  """
}

process CONVERT_LSD_NEXUS_TO_NEWICK {
  input:
  tuple val(size), val(rep), path(lsd_nexus)
  output:
  tuple val(size), val(rep), path("lsd_newick.nxs")

  """
  python $baseDir/scripts/helper.py 2 ${lsd_nexus} lsd_newick.nxs
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

  PREPARE_PHYLOTORCH(data.combine(use_bito))

  RUN_PHYLOTORCH(data.map { it.take(4) }.join(PREPARE_PHYLOTORCH.out, by: [0, 1]))
}
