#!/bin/bash

set -e
# source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

/usr/bin/time nextflow  \
    -C ./configs/rhino.config \
    run main.nf \
    --results "batch-results" \
    -profile rhino \
    -with-report "batch-results"/nextflow_report.html \
    -with-trace  "batch-results"/trace.txt \
    -work-dir "batch-results/work/" \
    -resume
