#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

/usr/bin/time nextflow  \
    -C ./configs/rhino.config \
    run main.nf \
    --results "$(date -I)-quoll-results" \
    -profile standard \
    -with-report ./q_output/nextflow_report.html \
    -with-trace ./q_output/trace.txt \
    -work-dir ./q_output/work/ \
    -resume
