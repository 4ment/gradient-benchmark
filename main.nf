#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.iterations = '5000'

flu_dataset = Channel.of(20, 50, 100, 200 , 500, 750, 1000, 1250, 1500, 2000)
subst_model = Channel.of("JC69", "GTR")

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
    input:
        path seq_file
        path tree_file
        path phylostan_stan
        path phylostan_kl
    output:
        path 'phylostan', emit: phylostan_out 
        
    """
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
                  --eta 0.00000001 \
                  --elbo_samples 1 \
                  --samples 1
    """
}

process PREPARE_PHYSHER {

    input:
        val size
        path seq_file
        path lsd_newick
    output:
        path "physher.json", emit: physher_json
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

    input:
        path json_file
        path seq_file
        path tree_file
    """
    physher $json_file
    """
}

process PREPARE_PHYLOTORCH {

    input:
        val size
        path seq_file
        path lsd_newick
    output:
        path "phylotorch.json", emit: phylotorch_json
    """
    python $baseDir/scripts/helper.py 3 \
                                      $seq_file \
                                      $lsd_newick \
                                      $dates_dir/H3N2_HA_2011_2013_${size}_0.lsd_dates.txt \
                                      $phylotorch_jc69_template phylotorch.json \
                                      ${params.iterations}
    """
}

process RUN_PHYLOTORCH {

    input:
        path json_file
        path seq_file
        path tree_file
    """
    phylotorch $json_file
    """
}

process CREATE_SUB_FILES {

    input:
        val size
    output:
        path "H3N2_HA_2011_2013_${size}_0.nwk", emit: subtree_dates_file
        path "H3N2_HA_2011_2013_${size}_0.fasta", emit: seq_file
        path "H3N2_HA_2011_2013_${size}_0.lsd_dates.txt", emit: new_dates_file
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
        path subtree_dates_file
        path new_dates_file
    output:
        path "lsd.out.date.nexus", emit: lsd_tree_dated // branch=time
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
        path lsd_nexus
    output:
        path 'lsd_newick.nxs', emit: lsd_newick

    """
    python $baseDir/scripts/helper.py 2 ${lsd_nexus} lsd_newick.nxs
    """
}

workflow {
    flu_dataset | CREATE_SUB_FILES
    
    RUN_LSD(CREATE_SUB_FILES.out[0], CREATE_SUB_FILES.out[2])
    
    CONVERT_LSD_NEXUS_TO_NEWICK(RUN_LSD.out.lsd_tree_dated)
    
    COMPILE_PHYLOSTAN()
    RUN_PHYLOSTAN(
       CREATE_SUB_FILES.out.seq_file,
       CONVERT_LSD_NEXUS_TO_NEWICK.out.lsd_newick,
       COMPILE_PHYLOSTAN.out[0],
       COMPILE_PHYLOSTAN.out[1])
    
    PREPARE_PHYSHER(
           flu_dataset,
           CREATE_SUB_FILES.out.seq_file,
           CONVERT_LSD_NEXUS_TO_NEWICK.out.lsd_newick)
    RUN_PHYSHER(
       PREPARE_PHYSHER.out,
       CREATE_SUB_FILES.out.seq_file,
       CONVERT_LSD_NEXUS_TO_NEWICK.out.lsd_newick)
    
    PREPARE_PHYLOTORCH(
           flu_dataset,
           CREATE_SUB_FILES.out.seq_file,
           CONVERT_LSD_NEXUS_TO_NEWICK.out.lsd_newick)
    RUN_PHYLOTORCH(
       PREPARE_PHYLOTORCH.out,
       CREATE_SUB_FILES.out.seq_file,
       CONVERT_LSD_NEXUS_TO_NEWICK.out.lsd_newick)
}