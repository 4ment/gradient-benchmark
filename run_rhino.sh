#!/bin/bash

set -e
# source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

/usr/bin/time nextflow  \
    -C ./configs/rhino.config \
    run main.nf \
    --results "batch-results-$(date -I)" \
    -profile rhino \
    -with-report "batch-results-$(date -I)"/nextflow_report.html \
    -with-trace  "batch-results-$(date -I)"/trace.txt \
    -work-dir "batch-results-$(date -I)/work/"
