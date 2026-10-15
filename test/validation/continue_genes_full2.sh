#!/usr/bin/env bash
set -euo pipefail
cd /shared/projects/metainvert/exogap
source /etc/profile.d/modules.sh
module load nextflow/25.10.2 singularity
export NXF_SINGULARITY_CACHEDIR=/shared/ifbstor1/projects/metainvert/exogap/.singularity
export APPTAINER_TMPDIR=/shared/projects/metainvert/exogap/.apptainer_tmp
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
export TMPDIR="$APPTAINER_TMPDIR"

# Wait for this earlier validation's pending tasks to finish storing their outputs.
prior_validation_pid="${1:-}"
if [[ -n "$prior_validation_pid" ]]; then
    while kill -0 "$prior_validation_pid" 2>/dev/null; do
        sleep 10
    done
fi

# Earlier tasks captured emapper's database warning instead of its package version.
# Repair only their version metadata; analysis outputs remain intact.
TASK_EGGNOG_VERSION=$(singularity exec .singularity/quay.io-biocontainers-eggnog-mapper-2.1.12--pyhdfd78af_1.img python -c "from importlib.metadata import version; print(version('eggnog-mapper'))")
export TASK_EGGNOG_VERSION
python - <<'PYVERSIONS'
import os
from pathlib import Path
root = Path('test/test_full2/cache/genes_annotation/functional_annotation/eggnogmapper')
for path in root.glob('*/*/versions.yml'):
    text = path.read_text()
    if 'error retrieving eggnog-mapper DB data' in text:
        path.with_suffix('.yml.invalid').write_text(text)
        path.write_text(text.splitlines()[0] + '\n    eggnog-mapper: ' + os.environ['TASK_EGGNOG_VERSION'] + '\n')
PYVERSIONS

set +e
nextflow -log test/validation/genes_full2_fixed.log run main.nf \
    -profile singularity,test_full2 -c conf/ifb.config \
    -c test/validation/slurm_runtime.config -ansi-log false
validation_status=$?
printf '%s\n' "$validation_status" > test/validation/genes_full2_fixed.exit
if [[ "$validation_status" -eq 0 ]]; then
    # Verify persistent output reuse without -resume.
    nextflow -log test/validation/genes_full2_store.log run main.nf \
        -profile singularity,test_full2 -c conf/ifb.config \
        -c test/validation/slurm_runtime.config -ansi-log false
    validation_status=$?
    printf '%s\n' "$validation_status" > test/validation/genes_full2_store.exit
fi
exit "$validation_status"
