#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reuse = false
params.results = "results"
params.enable_beast = false

include { treetime_validation } from "./modules/treetime_validation.nf" addParams(base: "$baseDir/treetime_validation")
include { micro } from "./modules/micro.nf"
include { macro_flu } from "./modules/macro_flu.nf"

dataset = "${baseDir}/treetime_validation/flu_H3N2/subtree_samples/dataset"
subtrees_alignment = "$baseDir/treetime_validation/resources/flu_H3N2/H3N2_HA_2011_2013.fasta"


process RUN_LSD {
  label 'ultrafast'

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
       -s 1701 \
       -l -1
  """
}

process CONVERT_LSD_NEXUS_TO_NEWICK {
  label 'ultrafast'

  input:
  tuple val(size), val(rep), path(lsd_nexus)
  output:
  tuple val(size), val(rep), path("lsd_newick.nxs")

  """
  helper.py 2 ${lsd_nexus} lsd_newick.nxs
  """
}

def group_per_size_rep(newick_ch, create_sub_ch) {
  newick_ch.join(
          create_sub_ch.map {
            v ->
              def size = v[0]
              def rep = v[1]
              def fasta = v[3]
              def dates = v[4]
              tuple(size, rep, fasta, dates)
          }, by: [0, 1])
}

process CREATE_SUB_FILES {
  label 'ultrafast'

  input:
  tuple val(size), val(rep), path(lsd_dates), path(newick_file)
  output:
  tuple val(size),
          val(rep),
          path("H3N2_HA_2011_2013_${size}_${rep}.new.nwk"),
          path("H3N2_HA_2011_2013_${size}_${rep}.fasta"),
          path("H3N2_HA_2011_2013_${size}_${rep}.lsd_dates.new.txt")
  """
  helper.py 0 \
            $subtrees_alignment \
            $lsd_dates \
            $newick_file \
            H3N2_HA_2011_2013_${size}_${rep}.new.nwk \
            H3N2_HA_2011_2013_${size}_${rep}.fasta \
            H3N2_HA_2011_2013_${size}_${rep}.lsd_dates.new.txt
  """
}


workflow {
  if (params.reuse) {
    subsets_ch = Channel.of(20, 50, 100, 200, 500, 750, 1000, 1250, 1500, 2000)
    replicates_ch = Channel.of(0..5)
    ch = subsets_ch.combine(replicates_ch)
    tt_ch = ch.map {
      tuple(it[0], it[1],
              file("${dataset}/LSD_out/H3N2_HA_2011_2013_${it[0]}_${it[1]}.lsd_dates.txt"),
              file("${dataset}/subtrees/H3N2_HA_2011_2013_${it[0]}_${it[1]}.nwk"))
    }
  } else {
    treetime_validation()
    tt_ch = treetime_validation.out
  }

  CREATE_SUB_FILES(tt_ch)

  RUN_LSD(CREATE_SUB_FILES.out)

  CONVERT_LSD_NEXUS_TO_NEWICK(RUN_LSD.out[0])

  data = group_per_size_rep(
          CONVERT_LSD_NEXUS_TO_NEWICK.out,
          CREATE_SUB_FILES.out)

  macro_flu(data)

  micro(data.map { it.take(4) })
}