#!/bin/bash
#SBATCH --job-name=annoglyph-myriapoda
#SBATCH --output=logs/myriapoda-%j.out
#SBATCH --error=logs/myriapoda-%j.err
#SBATCH --partition=long
#SBATCH --account=myriaphylo
#SBATCH --exclude=cpu-node-127,cpu-node-128,cpu-node-129,cpu-node-130
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00

set -euo pipefail
source /etc/profile.d/modules.sh
module load nextflow/25.10.2 singularity
ANNOGLYPH_PROJECT_DIR="${ANNOGLYPH_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$ANNOGLYPH_PROJECT_DIR"
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$ANNOGLYPH_PROJECT_DIR/.singularity}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$ANNOGLYPH_PROJECT_DIR/.apptainer_tmp}"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_TMPDIR"

test -s cache_save/cache_reuse_v1/manifest.json
test -s samplesheet_myriapoda_cache_save.csv
nextflow -log "logs/myriapoda-${SLURM_JOB_ID:-local}.nextflow.log" run main.nf -profile singularity -c conf/ifb.config \
    -c nextflow_myriapoda.config "$@"
