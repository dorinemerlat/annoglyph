#!/usr/bin/env bash
# Extend the published Python image with the statistical dependency used by reports.
set -euo pipefail
command -v singularity >/dev/null
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
base_image="${1:-$project_dir/.singularity/dorinemerlat-exogap_python-v1.5.img}"
output_image="${2:-$project_dir/.singularity/annoglyph_python-v1.6.sif}"
if [[ -e "$output_image" ]]; then
    echo "Output image already exists: $output_image" >&2
    exit 1
fi
build_parent="${APPTAINER_TMPDIR:-$project_dir/.apptainer_tmp}"
mkdir -p "$build_parent" "$(dirname "$output_image")"
export APPTAINER_TMPDIR="$build_parent"
export SINGULARITY_TMPDIR="$build_parent"
export TMPDIR="$build_parent"
if [[ ! -f "$base_image" ]]; then
    singularity pull "$base_image" docker://dorinemerlat/exogap_python:v1.5
fi
build_dir="$(mktemp -d "$build_parent/annoglyph-python.XXXXXXXX")"
trap 'rm -rf "$build_dir"' EXIT
singularity build --sandbox "$build_dir/rootfs" "$base_image"
mkdir -p "$build_dir/rootfs/shared"/{bank,projects,software,ifbstor1,home}
singularity exec --writable --pwd /tmp "$build_dir/rootfs" python -m pip install --no-cache-dir scipy==1.13.1
singularity exec --pwd /tmp "$build_dir/rootfs" python -c 'import Bio, ete3, gffpandas, gffutils, matplotlib, numpy, pandas, requests, scipy, seaborn; print("AnnoGlyph Python dependencies OK")'
singularity build --mksquashfs-args "-processors 4 -mem 2G" "$build_dir/annoglyph_python.sif" "$build_dir/rootfs"
mv "$build_dir/annoglyph_python.sif" "$output_image"
