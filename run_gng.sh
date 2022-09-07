#!/bin/bash

set -e
# source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

/usr/bin/time nextflow  \
    -C ./configs/rhino.config \
    run main.nf \
    --results "$(date -I)-gng-results" \
    -profile standard \
    -with-report "$(date -I)-gng-results"/nextflow_report.html \
    -with-trace  "$(date -I)-gng-results"/trace.txt \
    -work-dir '/fh/scratch/delete30/matsen_e/mathieu/temp/gng_work/'
