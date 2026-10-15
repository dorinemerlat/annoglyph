#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

HANSENIASPORA_DIR="${SCRIPT_DIR}/test_hanseniaspora"
PIERIS_DIR="${SCRIPT_DIR}/test_pieris"
UTRICULARIA_DIR="${SCRIPT_DIR}/test_utricularia"

HANSENIASPORA_ACCESSIONS=(
    "GCA_050947685.1" # Hanseniaspora vineae
    "GCA_050947655.1" # Hanseniaspora osmophila
    "GCA_050947605.1" # Hanseniaspora guilliermondii
    "GCA_050947625.1" # Hanseniaspora opuntiae
    "GCA_040256855.1" # Hanseniaspora menglaensis
)

PIERIS_ACCESSIONS=(
    "GCF_905147105.1" # Pieris brassicae
    "GCF_905475465.1" # Pieris napi
    "GCF_905147795.1" # Pieris rapae
)

UTRICULARIA_ACCESSIONS=(
    "GCA_982292365.1" # Utricularia gibba
    "GCA_056322005.1" # Utricularia aurea
    "GCA_009725065.1" # Utricularia reniformis
)

download_set() {
    local output_dir="$1"
    shift
    local accessions=("$@")

    rm -rf "${output_dir}"
    mkdir -p "${output_dir}"

    echo "Downloading ${#accessions[@]} assemblies to ${output_dir}..."

    datasets download genome accession "${accessions[@]}" \
        --include genome,gff3,protein \
        --filename "${output_dir}/genomes.zip"

    unzip -q "${output_dir}/genomes.zip" -d "${output_dir}"
    rm "${output_dir}/genomes.zip"

    for accession in "${accessions[@]}"; do
        mv "${output_dir}/ncbi_dataset/data/${accession}" "${output_dir}/"
    done

    rm -rf "${output_dir}/ncbi_dataset"
    rm -f "${output_dir}/README.md" "${output_dir}/md5sum.txt"
}

download_set "${HANSENIASPORA_DIR}" "${HANSENIASPORA_ACCESSIONS[@]}"
download_set "${PIERIS_DIR}" "${PIERIS_ACCESSIONS[@]}"
download_set "${UTRICULARIA_DIR}" "${UTRICULARIA_ACCESSIONS[@]}"

echo
echo "Test datasets downloaded:"
echo "  ${HANSENIASPORA_DIR}"
echo "  ${PIERIS_DIR}"
echo "  ${UTRICULARIA_DIR}"

# Generate samplesheets and small clade-specific protein evidence for Nextflow.
python3 "${SCRIPT_DIR}/prepare_test_sets.py"
