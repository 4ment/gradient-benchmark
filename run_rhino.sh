#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

/usr/bin/time nextflow  \
    -C ./configs/rhino.config \
    run main.nf \
    -profile rhino \
    -with-report ./r_output/nextflow_report.html \
    -work-dir ./r_output/work/ \
    -resume
