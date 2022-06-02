
nextflow.enable.dsl = 2

process COMBIME_TIME_LOG {
  label 'ultrafast'

  publishDir "$params.results/$tag/", mode: 'copy'

  input:
    path files
    val(tag)
  output:
    path("${tag}.csv")
  """
  time-parser.py ${tag}.csv ${files}
  """
}

process CONVERT_LSD_NEXUS_TO_NEWICK {
  label 'ultrafast'

  input:
    tuple val(size), val(rep), path(lsd_nexus), val(rate)
  output:
    tuple val(size), val(rep), path('lsd_newick.nxs'), val(rate)
  """
  helper.py convert --input ${lsd_nexus} --output lsd_newick.nxs
  """
}

process CREATE_SUB_FILES {
  label 'ultrafast'

  input:
    tuple val(size), val(rep), path(lsd_dates), path(newick_file), path(fasta_file)
  output:
    tuple val(size), val(rep),
          path("H3N2_HA_2011_2013_${size}_${rep}.new.nwk"),
          path("H3N2_HA_2011_2013_${size}_${rep}.fasta"),
          path("H3N2_HA_2011_2013_${size}_${rep}.lsd_dates.new.txt")
  """
  helper.py subfiles \
            --input $fasta_file \
            --dates $lsd_dates \
            --tree $newick_file \
            --out_tree H3N2_HA_2011_2013_${size}_${rep}.new.nwk \
            --out_fasta H3N2_HA_2011_2013_${size}_${rep}.fasta \
            --out_dates H3N2_HA_2011_2013_${size}_${rep}.lsd_dates.new.txt
  """
}
