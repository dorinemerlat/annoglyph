#!/bin/bash
#SBATCH --job-name=exogap
#SBATCH --output=exogap-%j.out
#SBATCH --error=exogap-%j.err
#SBATCH --mail-type ALL
#SBATCH --mail-user dorine.merlat@etu.unistra.fr
#SBATCH --partition fast
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks (Nextflow runs as a single task)
#SBATCH --account=orthoexplorer
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00

set -euo pipefail
source /etc/profile.d/modules.sh
module load nextflow/25.10.2 singularity
EXOGAP_PROJECT_DIR="${EXOGAP_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$EXOGAP_PROJECT_DIR/.singularity}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$EXOGAP_PROJECT_DIR/.apptainer_tmp}"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
export TMPDIR="$APPTAINER_TMPDIR"

cd "$EXOGAP_PROJECT_DIR"
mkdir -p "$APPTAINER_TMPDIR"

# Usage: sbatch nextflow.sbatch test_hanseniaspora_minimal [additional Nextflow options]
test_profile="${1:-test_hanseniaspora_minimal}"
if (( $# > 0 )); then shift; fi
case "$test_profile" in
    test_hanseniaspora|test_pieris|test_utricularia|test_hanseniaspora_minimal|test_pieris_minimal|test_utricularia_minimal|test_hanseniaspora_full|test_pieris_full|test_utricularia_full|test|test_full1|test_full2) ;;
    *) echo "Unknown test profile: $test_profile" >&2; exit 2 ;;
esac

nextflow run main.nf -profile "singularity,$test_profile" -c conf/ifb.config \
    -c test/validation/slurm_runtime.config "$@"
