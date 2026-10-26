#!/usr/bin/env bash
#SBATCH --job-name=annoglyph-repeat-assets
#SBATCH --partition=long
#SBATCH --account=myriaphylo
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=logs/finalize-repeat-assets-%j.out
#SBATCH --error=logs/finalize-repeat-assets-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
python bin/publish_cached_repeat_library.py
