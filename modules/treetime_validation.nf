#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.subtrees = Channel.of(20, 50, 100, 200, 500, 750, 1000)//, 1250, 1500, 2000)
params.subtrees_replicates = Channel.of(0..9)

if(params.enable_beast){
  beast_res = "--beast_file beast_res.csv"
}else{
  beast_res = ""
}

params.base = "$baseDir"

alignment_file = "$params.base/resources/flu_H3N2/H3N2_HA_2011_2013.fasta"
tree_file = "$params.base/resources/flu_H3N2/H3N2_HA_2011_2013.nwk"
beast_template = "$params.base/resources/beast/template_bedford_et_al_2015.xml"
treetime_flu_H3N2 = "$params.base/flu_H3N2/subtree_samples"

process TREETIME_VALIDATION_SUBTREES {
  label 'treetime'

  publishDir "${treetime_flu_H3N2}", mode: 'copy'

  input:
  tuple val(size), val(rep)
  output:
  tuple val(size),
          val(rep),
          path("dataset/LSD_out/H3N2_HA_2011_2013_${size}_${rep}.lsd_dates.txt"),
          path("dataset/subtrees/H3N2_HA_2011_2013_${size}_${rep}.nwk")
  tuple path("treetime_res.${size}.${rep}.csv"),
          path("lsd_res.${size}.${rep}.csv")
  path("beast_res.${size}.${rep}.csv") optional true
  script:
  if(params.enable_beast)
    beast_res = "--beast_file beast_res.${size}.${rep}.csv"
  else
    beast_res = ""
  
  """
  python2.7 $params.base/generate_flu_subtrees_dataset_run.py --size $size \
                                                           --out_dir dataset \
                                                           --suffix $rep \
                                                           --treetime_file treetime_res.${size}.${rep}.csv \
                                                           --lsd_file lsd_res.${size}.${rep}.csv \
                                                           --aln_file ${alignment_file} \
                                                           --tree_file ${tree_file} \
                                                           --template_file ${beast_template} \
                                                           $beast_res
  """

}

process COMBIME_CSV {
  label 'ultrafast'

  publishDir "${treetime_flu_H3N2}", mode: 'copy'

  input:
  path files
  output:
  path("lsd_res.csv")
  path("treetime_res.csv")
  path("beast_res.csv") optional true
  """
  head -n1 lsd_res.20.0.csv > lsd_res.csv
  tail -q -n+2 lsd_res*[0-9].csv >> lsd_res.csv
  head -n1 treetime_res.20.0.csv > treetime_res.csv
  tail -q -n+2 treetime_res*[0-9].csv >> treetime_res.csv
  if [ -f "beast_res.20.0.csv" ]; then
    head -n1 beast_res.20.0.csv > beast_res.csv
    tail -q -n+2 beast_res*[0-9].csv >> beast_res.csv
  fi
  """
}

workflow treetime_validation {
  main:
  params.subtrees.combine(params.subtrees_replicates) | TREETIME_VALIDATION_SUBTREES
  
  ch_csv = TREETIME_VALIDATION_SUBTREES.out[1].collect()

  if (params.enable_beast)
    ch_csv = ch_csv.mix(TREETIME_VALIDATION_SUBTREES.out[2].collect())

  COMBIME_CSV(ch_csv.collect())

  emit:
  TREETIME_VALIDATION_SUBTREES.out[0]
}
