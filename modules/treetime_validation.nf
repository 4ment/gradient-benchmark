#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.subtrees = Channel.of(20, 50, 100, 200, 500, 750, 1000, 1250, 1500, 2000)
params.subtrees_replicates = Channel.of(0..9)

params.base = "$baseDir"

alignment_file = "$params.base/resources/flu_H3N2/H3N2_HA_2011_2013.fasta"
tree_file = "$params.base/resources/flu_H3N2/H3N2_HA_2011_2013.nwk"
beast_template = "$params.base/resources/beast/template_bedford_et_al_2015.xml"
treetime_flu_H3N2 = "$params.base/flu_H3N2/subtree_samples"

process TREETIME_VALIDATION_SUBTREES {
  publishDir "${treetime_flu_H3N2}", mode: 'copy'

  input:
  tuple val(size), val(rep)
  output:
  tuple val(size),
          val(rep),
          path("dataset/LSD_out/H3N2_HA_2011_2013_${size}_${rep}.lsd_dates.txt"),
          path("dataset/subtrees/H3N2_HA_2011_2013_${size}_${rep}.nwk")
  path("treetime_res.csv")
  path("lsd_res.csv")
  path("beast_res.csv")
  """
  source activate treetime
  python $params.base/generate_flu_subtrees_dataset_run.py $size \
                                                           dataset \
                                                           $rep \
                                                           treetime_res.csv \
                                                           lsd_res.csv \
                                                           beast_res.csv \
                                                           ${alignment_file} \
                                                           ${tree_file} \
                                                           ${beast_template}
  """

}

workflow treetime_validation {
  main:
  params.subtrees.combine(params.subtrees_replicates) | TREETIME_VALIDATION_SUBTREES
  emit:
  TREETIME_VALIDATION_SUBTREES.out[0]
}