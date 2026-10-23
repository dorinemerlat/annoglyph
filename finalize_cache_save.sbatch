#!/usr/bin/env bash
#SBATCH --job-name=annoglyph-finalize
#SBATCH --partition=long
#SBATCH --account=myriaphylo
#SBATCH --exclude=cpu-node-127,cpu-node-128,cpu-node-129,cpu-node-130
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/finalize-%j.out
#SBATCH --error=logs/finalize-%j.err
set -euo pipefail
source /etc/profile.d/modules.sh
module load nextflow/25.10.2 singularity
cd "${SLURM_SUBMIT_DIR:?}"
export NXF_SINGULARITY_CACHEDIR="$PWD/.singularity"
export APPTAINER_TMPDIR="$PWD/.apptainer_tmp"
mkdir -p "$APPTAINER_TMPDIR"
nextflow -log "logs/finalize-${SLURM_JOB_ID}.nextflow.log" run finalize_cache_save.nf \
    -profile singularity -c conf/ifb.config -c nextflow_myriapoda.config \
    -c conf/finalize_cache_save.config "$@"
