#!/bin/bash

module load system/singularity
module load devel/miniforge
conda activate nextflow

mkdir -p ./scratch
export SCRATCH="$(pwd)/scratch"

# Replace this with the name of your samplesheet
samplesheet="./design.csv"
rnaseq_version="3.20.0"
genome="hg38"
profile="${HOME}/helix.config"

# PARAMETERS SECTION
# ==========================================================


# Replace this with the name of your samplesheet
samplesheet="./samplesheet.csv"

aligner="star_rsem"
outdir="results"

# NF-CORE RUN
# ===========================================================
 
#~/matrix/matrix.sh "Start nf-core job $SLURM_JOB_ID"

nextflow run nf-core/rnaseq \
    -c "${profile}" \
    -r "${rnaseq_version}" \
    -bg \
    -resume \
    --input "${samplesheet}" \
    --genome "${genome}" \
    --aligner "${aligner}" \
    --save-align-intermeds \
    --outdir "${outdir}" > nextflow.run.2.out 2> nextflow.run.2.err

