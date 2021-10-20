#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.results = "results"

params.subtrees_alignment = "$baseDir/treetime_validation/resources/flu_H3N2/H3N2_HA_2011_2013.fasta"

include { treetime_validation } from "./modules/treetime_validation.nf" addParams(base: "$baseDir/treetime_validation")
include { micro } from "./modules/micro.nf"
include { macro_flu } from "./modules/macro_flu.nf"

process RUN_LSD {
  input:
  tuple val(size),
          val(rep),
          path(subtree_dates_file),
          path(fasta),
          path(new_dates_file)
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

def group_per_size(newick_ch, create_sub_ch) {
  newick_ch
          .phase(create_sub_ch)
          .map { left, right ->
            def size = left[0]
            def rep = left[1]
            def tree = left[2]
            def fasta = right[3]
            tuple(size, rep, tree, fasta)
          }
}

def group_per_size2(newick_ch, create_sub_ch) {
  newick_ch
          .phase(create_sub_ch)
          .map { left, right ->
            def size = left[0]
            def rep = left[1]
            def tree = left[2]
            def dates = right[4]
            def fasta = right[3]
            tuple(size, rep, tree, fasta, dates)
          }
}


process CREATE_SUB_FILES {

  input:
  tuple val(size), val(rep), path(lsd_dates), path(newick_file)
  output:
  tuple val(size),
          val(rep),
          path("H3N2_HA_2011_2013_${size}_${rep}.nwk"),
          path("H3N2_HA_2011_2013_${size}_${rep}.fasta"),
          path("H3N2_HA_2011_2013_${size}_${rep}.lsd_dates.txt")
  """
  python $baseDir/scripts/helper.py 0 \
                                    $params.subtrees_alignment \
                                    $lsd_dates \
                                    $newick_file \
                                    H3N2_HA_2011_2013_${size}_${rep}.nwk \
                                    H3N2_HA_2011_2013_${size}_${rep}.fasta \
                                    H3N2_HA_2011_2013_${size}_${rep}.lsd_dates.txt
  """
}


workflow {
  treetime_validation()

  CREATE_SUB_FILES(treetime_validation.out)

  RUN_LSD(CREATE_SUB_FILES.out)

  CONVERT_LSD_NEXUS_TO_NEWICK(RUN_LSD.out[0])

  macro_flu(group_per_size2(
          CONVERT_LSD_NEXUS_TO_NEWICK.out,
          CREATE_SUB_FILES.out))

  micro(group_per_size(
          CONVERT_LSD_NEXUS_TO_NEWICK.out,
          CREATE_SUB_FILES.out))
}