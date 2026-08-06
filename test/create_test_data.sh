#!/usr/bin/env bash
set -euo pipefail

module load ncbi-datasets-cli/18.14.0
module load seqkit/2.13.0

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

TEST_FULL1_DIR="${SCRIPT_DIR}/test_full1"
TEST_FULL2_DIR="${SCRIPT_DIR}/test_full2"
TEST_MINIMAL_DIR="${SCRIPT_DIR}/test_minimal"
DOWNLOAD_DIR="${SCRIPT_DIR}/.test_data_download"

FULL1_ACCESSIONS=(
    "GCA_001747055.1"
    "GCA_037102615.1"
    "GCA_050947715.1"
)

FULL2_ACCESSIONS=(
    "GCF_000146045.2"
    "GCF_002079055.1"
)

MINIMAL_ACCESSION="GCA_050947715.1"
MINIMAL_SEQUENCES=(
    "CM116924.1"
    "CM116925.1"
    "CM116926.1"
)

ALL_ACCESSIONS=(
    "${FULL1_ACCESSIONS[@]}"
    "${FULL2_ACCESSIONS[@]}"
)

rm -rf \
    "${DOWNLOAD_DIR}" \
    "${TEST_FULL1_DIR}" \
    "${TEST_FULL2_DIR}" \
    "${TEST_MINIMAL_DIR}"

mkdir -p \
    "${DOWNLOAD_DIR}" \
    "${TEST_FULL1_DIR}/genomes" \
    "${TEST_FULL2_DIR}/genomes" \
    "${TEST_MINIMAL_DIR}/genomes"

# Download genome assemblies and reference annotations.
datasets download genome accession \
    "${ALL_ACCESSIONS[@]}" \
    --include genome,gff3,protein,cds \
    --filename "${DOWNLOAD_DIR}/genomes.zip"

unzip -q "${DOWNLOAD_DIR}/genomes.zip" -d "${DOWNLOAD_DIR}"

copy_assembly_files() {
    local accession=$1
    local test_dir=$2
    local source_dir="${DOWNLOAD_DIR}/ncbi_dataset/data/${accession}"
    local target_dir="${test_dir}/genomes/${accession}"
    local genome

    mkdir -p "${target_dir}"

    genome=$(find "${source_dir}" -maxdepth 1 -type f -name "${accession}_*_genomic.fna" -print -quit)

    if [[ -z "${genome}" ]]; then
        echo "Genome FASTA not found for ${accession}" >&2
        exit 1
    fi

    for required_file in genomic.gff protein.faa cds_from_genomic.fna; do
        if [[ ! -f "${source_dir}/${required_file}" ]]; then
            echo "Required file not found for ${accession}: ${required_file}" >&2
            exit 1
        fi
    done

    cp "${genome}" "${target_dir}/genome.fa"
    cp "${source_dir}/genomic.gff" "${target_dir}/annotation.gff"
    cp "${source_dir}/protein.faa" "${target_dir}/proteins.faa"
    cp "${source_dir}/cds_from_genomic.fna" "${target_dir}/cds.fna"
}

# Create the first full test dataset with Hanseniaspora uvarum assemblies.
for accession in "${FULL1_ACCESSIONS[@]}"; do
    copy_assembly_files "${accession}" "${TEST_FULL1_DIR}"
done

# Create the second full test dataset with Saccharomyces species.
for accession in "${FULL2_ACCESSIONS[@]}"; do
    copy_assembly_files "${accession}" "${TEST_FULL2_DIR}"
done

# Extract three chromosomes for the minimal test dataset.
minimal_source="${TEST_FULL1_DIR}/genomes/${MINIMAL_ACCESSION}"
minimal_fasta="${TEST_MINIMAL_DIR}/genomes/hanseniaspora_uvarum_chr1-3.fa"
minimal_gff="${TEST_MINIMAL_DIR}/genomes/hanseniaspora_uvarum_chr1-3.gff"
sequence_list="${DOWNLOAD_DIR}/minimal_sequences.txt"

printf "%s\n" "${MINIMAL_SEQUENCES[@]}" > "${sequence_list}"

seqkit grep \
    --pattern-file "${sequence_list}" \
    "${minimal_source}/genome.fa" \
    > "${minimal_fasta}"

awk -v sequence_file="${sequence_list}" '
BEGIN {
    FS = OFS = "\t"

    while ((getline sequence_id < sequence_file) > 0) {
        selected[sequence_id] = 1
    }

    close(sequence_file)
}

/^##gff-version/ {
    print
    next
}

/^##sequence-region/ {
    if ($2 in selected) {
        print
    }
    next
}

/^#/ {
    next
}

$1 in selected {
    print
}
' "${minimal_source}/annotation.gff" > "${minimal_gff}"

# Validate that all requested chromosomes were extracted.
extracted_sequence_count=$(grep -c '^>' "${minimal_fasta}")

if [[ "${extracted_sequence_count}" -ne "${#MINIMAL_SEQUENCES[@]}" ]]; then
    echo "Expected ${#MINIMAL_SEQUENCES[@]} sequences but extracted ${extracted_sequence_count}." >&2
    echo "Requested sequences: ${MINIMAL_SEQUENCES[*]}" >&2
    exit 1
fi

# Validate the generated datasets.
echo "Full test dataset 1:"
seqkit stats "${TEST_FULL1_DIR}"/genomes/*/genome.fa

echo
echo "Full test dataset 2:"
seqkit stats "${TEST_FULL2_DIR}"/genomes/*/genome.fa

echo
echo "Minimal test dataset:"
seqkit stats "${minimal_fasta}"

minimal_gene_count=$(awk -F'\t' '$3 == "gene" {count++} END {print count + 0}' "${minimal_gff}")

echo
echo "Minimal test genome contains ${extracted_sequence_count} sequences."
echo "Minimal test annotation contains ${minimal_gene_count} genes."

# Remove temporary NCBI download files.
rm -rf "${DOWNLOAD_DIR}"