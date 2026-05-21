#!/usr/bin/env bash

set -euo pipefail

MAP="${MAP:-/gpfs/chencao/qinminzhang/Nature/mine/code/configs/path.file_id_map.tsv}"
K="${K:-60}"

ESSENTIAL_DIR="${ESSENTIAL_DIR:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/run_all/outputs/perturbseq/cnmf_essential/trait_association/K562_essential_raw_singlecell_01/ProgramLevel}"
GENOMEWIDE_DIR="${GENOMEWIDE_DIR:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/run_all/outputs/perturbseq/cnmf_genomewide/trait_association/K562GW/ProgramLevel}"

if [ ! -f "${MAP}" ]; then
    echo "MAP not found: ${MAP}" >&2
    exit 1
fi

map_traits=$(( $(wc -l < "${MAP}") - 1 ))
expected_files_per_dataset=$((map_traits * 2))
essential_actual_files=$(find "${ESSENTIAL_DIR}" -maxdepth 1 -type f 2>/dev/null | wc -l || true)
genomewide_actual_files=$(find "${GENOMEWIDE_DIR}" -maxdepth 1 -type f 2>/dev/null | wc -l || true)

printf 'input_summary\n'
printf 'dataset\ttraits_in_map\texpected_files\tactual_files\n'
printf 'cnmf_essential\t%s\t%s\t%s\n' "${map_traits}" "${expected_files_per_dataset}" "${essential_actual_files}"
printf 'cnmf_genomewide\t%s\t%s\t%s\n' "${map_traits}" "${expected_files_per_dataset}" "${genomewide_actual_files}"
printf '\nmissing_files\n'
printf 'dataset\tid2\ttrait_file\tmissing_programs\tmissing_regulators\n'

total_missing=0
essential_missing=0
genomewide_missing=0
essential_missing_programs=0
essential_missing_regulators=0
genomewide_missing_programs=0
genomewide_missing_regulators=0
expected_pairs=0

find_existing_file() {
    local dir="$1"
    local prefix="$2"
    shift 2

    local candidate
    for candidate in "$@"; do
        if [ -f "${dir}/${prefix}${candidate}" ]; then
            printf '%s' "${candidate}"
            return 0
        fi
    done
    return 1
}

while IFS=$'\t' read -r id1 id2 path1 path2; do
    path2="${path2%$'\r'}"
    source_file="$(basename "${path2}")"
    trait_file="${source_file}.per_gene_estimates.tsv"
    legacy_trait_file="${source_file}"

    for item in \
        "cnmf_essential:${ESSENTIAL_DIR}" \
        "cnmf_genomewide:${GENOMEWIDE_DIR}"
    do
        dataset="${item%%:*}"
        dir="${item#*:}"
        expected_pairs=$((expected_pairs + 1))

        missing_programs=0
        missing_regulators=0
        if ! find_existing_file "${dir}" "programs_enrichment_K${K}_" "${trait_file}" "${legacy_trait_file}" >/dev/null; then
            missing_programs=1
        fi
        if ! find_existing_file "${dir}" "regulators_enrichment_K${K}_" "${trait_file}" "${legacy_trait_file}" >/dev/null; then
            missing_regulators=1
        fi

        if [ "${missing_programs}" -ne 0 ] || [ "${missing_regulators}" -ne 0 ]; then
            printf '%s\t%s\t%s\t%s\t%s\n' \
                "${dataset}" \
                "${id2}" \
                "${trait_file}" \
                "${missing_programs}" \
                "${missing_regulators}"

            total_missing=$((total_missing + 1))
            if [ "${dataset}" = "cnmf_essential" ]; then
                essential_missing=$((essential_missing + 1))
                essential_missing_programs=$((essential_missing_programs + missing_programs))
                essential_missing_regulators=$((essential_missing_regulators + missing_regulators))
            else
                genomewide_missing=$((genomewide_missing + 1))
                genomewide_missing_programs=$((genomewide_missing_programs + missing_programs))
                genomewide_missing_regulators=$((genomewide_missing_regulators + missing_regulators))
            fi
        fi
    done
done < <(tail -n +2 "${MAP}")

printf '\nsummary\n'
printf 'dataset\tmissing_traits\tmissing_programs\tmissing_regulators\n'
printf 'cnmf_essential\t%s\t%s\t%s\n' \
    "${essential_missing}" \
    "${essential_missing_programs}" \
    "${essential_missing_regulators}"
printf 'cnmf_genomewide\t%s\t%s\t%s\n' \
    "${genomewide_missing}" \
    "${genomewide_missing_programs}" \
    "${genomewide_missing_regulators}"
printf 'all\t%s\t%s\t%s\n' \
    "${total_missing}" \
    "$((essential_missing_programs + genomewide_missing_programs))" \
    "$((essential_missing_regulators + genomewide_missing_regulators))"
printf 'expected_dataset_trait_pairs\t%s\n' "${expected_pairs}"
