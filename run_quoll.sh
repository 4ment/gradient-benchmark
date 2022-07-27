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
    -with-report "$(date -I)-quoll-results"/nextflow_report.html \
    -with-trace  "$(date -I)-quoll-results"/trace.txt \
    -work-dir ./q_work/
