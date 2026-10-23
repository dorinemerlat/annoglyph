#!/usr/bin/env bash
#SBATCH --job-name=annoglyph-validate
#SBATCH --partition=long
#SBATCH --account=myriaphylo
#SBATCH --exclude=cpu-node-127,cpu-node-128,cpu-node-129,cpu-node-130
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=logs/finalize-validation-%j.out
#SBATCH --error=logs/finalize-validation-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
python bin/validate_cached_finalization.py
