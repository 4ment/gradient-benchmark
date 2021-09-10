#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.iterations = '5000'
params.results    = "results"

flu_dataset = Channel.of(20, 50, 100, 200 , 500, 750, 1000, 1250, 1500, 2000)
subst_model = Channel.of("JC69", "GTR")
use_bito = Channel.of(true, false)

flu_H3N2 = "$baseDir/flu_H3N2"
subdata_dir = "$flu_H3N2/results"
physher_jc69_template = "$flu_H3N2/physher-JC69.template"
phylotorch_jc69_template = "$flu_H3N2/phylotorch-JC69.template"

// treetime directory
treetime_base = "$baseDir/treetime_validation"
flu_H3N2_dataset = "${treetime_base}/flu_H3N2/subtree_samples/dataset"
dates_dir = "${flu_H3N2_dataset}/LSD_out"
subtrees_dir = "${flu_H3N2_dataset}/subtrees"

alignment_file = "$treetime_base/resources/flu_H3N2/H3N2_HA_2011_2013.fasta"
tree_file = "$treetime_base/resources/flu_H3N2/H3N2_HA_2011_2013.nwk"


process COMPILE_PHYLOSTAN {
    output:
        path "H3N2_HA_2011_2013.stan", emit: phylostan_stan
        path "H3N2_HA_2011_2013.pkl", emit: phylostan_pkl
    """
    phylostan build -s H3N2_HA_2011_2013.stan \
    -m JC69 \
    --heterochronous \
    --estimate_rate \
    --clock strict \
    -c constant \
    --compile
    """
}

process RUN_PHYLOSTAN {
    publishDir "$params.results/phylostan.${size}"

    input:
        tuple val(size), path(tree_file), path(seq_file)
        path phylostan_stan
        path phylostan_kl
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
                  --samples 1 > out.txt ; }2 > time.log
    """
}

process PREPARE_PHYSHER {

    input:
        tuple val(size), path(lsd_newick), path(seq_file)
    output:
        tuple val(size), path("physher.json")
    """
    python $baseDir/scripts/helper.py 1 \
                                      $seq_file \
                                      $lsd_newick  \
                                      $dates_dir/H3N2_HA_2011_2013_${size}_0.lsd_dates.txt \
                                      $physher_jc69_template physher.json \
                                      ${params.iterations}
    """
}

process RUN_PHYSHER {
    publishDir "$params.results/physher.${size}"

    input:
        tuple val(size), path(lsd_newick), path(seq_file), path(physher_json)
    output:
        path("out.txt")
        path("time.log")
    """
    { time physher $physher_json > out.txt ; } 2> time.log
    """
}

process PREPARE_PHYLOTORCH {

    input:
        tuple val(size), path(lsd_newick), path(seq_file), val(bito)
    output:
        tuple val(size), path("phylotorch.json"), val(bito)
    """
    python $baseDir/scripts/helper.py 3 \
                                      $seq_file \
                                      $lsd_newick \
                                      $dates_dir/H3N2_HA_2011_2013_${size}_0.lsd_dates.txt \
                                      $phylotorch_jc69_template phylotorch.json \
                                      ${params.iterations} \
                                      ${bito}
    """
}

process RUN_PHYLOTORCH {
    publishDir "$params.results/phylotorch.${bito}.${size}"

    input:
        tuple val(size), path(lsd_newick), path(seq_file), path(phylotorch_json), val(bito)
    output:
        path("out.txt")
        path("time.log")
    """
    { time phylotorch $phylotorch_json > out.txt ; } 2> time.log
    """
}

process CREATE_SUB_FILES {

    input:
        val size
    output:
        tuple val(size), \
          path("H3N2_HA_2011_2013_${size}_0.nwk"), \
          path("H3N2_HA_2011_2013_${size}_0.fasta"), \
          path("H3N2_HA_2011_2013_${size}_0.lsd_dates.txt")
    """
    python $baseDir/scripts/helper.py 0 \
                                      $alignment_file \
                                      $dates_dir/H3N2_HA_2011_2013_${size}_0.lsd_dates.txt \
                                      $subtrees_dir/H3N2_HA_2011_2013_${size}_0.nwk \
                                      H3N2_HA_2011_2013_${size}_0.nwk \
                                      H3N2_HA_2011_2013_${size}_0.fasta  \
                                      H3N2_HA_2011_2013_${size}_0.lsd_dates.txt
    """
}

process RUN_LSD {
    input:
        tuple val(size), path(subtree_dates_file), path(fasta), path(new_dates_file)
    output:
        tuple val(size), path("lsd.out.date.nexus")
        path "lsd.out.nexus", emit: lsd_tree // branch=subst not used
        path "lsd.out.nwk", emit: lsd_tree_newick // branch=subst
        path "lsd.out"
        
    """
    lsd2 -i ${subtree_dates_file} \
         -d ${new_dates_file} \
         -o lsd.out \
         -s 1701 \
         -c
    """
}

process CONVERT_LSD_NEXUS_TO_NEWICK {
    input:
        tuple val(size), path(lsd_nexus)
    output:
        tuple val(size), path("lsd_newick.nxs")

    """
    python $baseDir/scripts/helper.py 2 ${lsd_nexus} lsd_newick.nxs
    """
}

def group_per_size(newick_ch, create_sub_ch) {
  newick_ch
    .phase(create_sub_ch)
    .map{ left, right ->
      def size = left[0]
      def tree = left[1]
      def fasta = right[2]
      tuple(size, tree, fasta)
    }
}

workflow {
    flu_dataset | CREATE_SUB_FILES
    
    RUN_LSD(CREATE_SUB_FILES.out)
    
    CONVERT_LSD_NEXUS_TO_NEWICK(RUN_LSD.out[0])
    
    COMPILE_PHYLOSTAN()
    RUN_PHYLOSTAN(
       group_per_size(
            CONVERT_LSD_NEXUS_TO_NEWICK.out,
            CREATE_SUB_FILES.out),
       COMPILE_PHYLOSTAN.out[0],
       COMPILE_PHYLOSTAN.out[1] )
    
    PREPARE_PHYSHER(
           group_per_size(
                CONVERT_LSD_NEXUS_TO_NEWICK.out,
                CREATE_SUB_FILES.out) )
    RUN_PHYSHER(
       group_per_size(
            CONVERT_LSD_NEXUS_TO_NEWICK.out,
            CREATE_SUB_FILES.out).join(PREPARE_PHYSHER.out) )

    
    PREPARE_PHYLOTORCH(
           group_per_size(
                CONVERT_LSD_NEXUS_TO_NEWICK.out,
                CREATE_SUB_FILES.out).combine(use_bito) )

    RUN_PHYLOTORCH(
       group_per_size(
            CONVERT_LSD_NEXUS_TO_NEWICK.out,
            CREATE_SUB_FILES.out).join(PREPARE_PHYLOTORCH.out) )
}