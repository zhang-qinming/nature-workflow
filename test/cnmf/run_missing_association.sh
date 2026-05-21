#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"

MAP="${MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"
K="${K:-60}"
RANDOM_ITERATIONS="${RANDOM_ITERATIONS:-100000}"
TOP_N="${TOP_N:-100}"

RUN_OUTPUT_ROOT="${RUN_OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/run_all/outputs}"
CODE_OUTPUT_ROOT="${CODE_OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/Nature/mine/code/outputs}"
LOF_DIR="${LOF_DIR:-${RUN_OUTPUT_ROOT}/genebayes/posterior}"

GENE_MAP="${GENE_MAP:-${PROJECT_ROOT}/data/gencode_v41_gname_gid_ALL_sorted_onlyID}"
SHET_PATH="${SHET_PATH:-${PROJECT_ROOT}/data/shet_10bins.txt}"
R_SCRIPT="${R_SCRIPT:-${PROJECT_ROOT}/src/paper_pipeline/dataprep/perturbseq/r_scripts/cnmf_burden_program_regulators.R}"

ESSENTIAL_CELL="${ESSENTIAL_CELL:-K562_essential_raw_singlecell_01}"
GENOMEWIDE_CELL="${GENOMEWIDE_CELL:-K562GW}"

ESSENTIAL_ASSOC_DIR="${ESSENTIAL_ASSOC_DIR:-${RUN_OUTPUT_ROOT}/perturbseq/cnmf_essential/trait_association/${ESSENTIAL_CELL}/ProgramLevel}"
GENOMEWIDE_ASSOC_DIR="${GENOMEWIDE_ASSOC_DIR:-${RUN_OUTPUT_ROOT}/perturbseq/cnmf_genomewide/trait_association/${GENOMEWIDE_CELL}/ProgramLevel}"

BATCH_ROOT="${BATCH_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/run_all/scripts/cnmf_association_rerun}"
JOBS_DIR="${JOBS_DIR:-${BATCH_ROOT}/jobs}"
LOGS_DIR="${LOGS_DIR:-${BATCH_ROOT}/logs}"
STATUS_DIR="${STATUS_DIR:-${BATCH_ROOT}/status}"
MANIFEST_PATH="${MANIFEST_PATH:-${BATCH_ROOT}/manifest.tsv}"

SBATCH_PARTITION="${SBATCH_PARTITION:-fat}"
SBATCH_CPUS="${SBATCH_CPUS:-10}"
SBATCH_MEM="${SBATCH_MEM:-160G}"
SBATCH_TIME="${SBATCH_TIME:-24:00:00}"
JOB_NAME="${JOB_NAME:-cnmf_assoc}"
DRY_RUN="${DRY_RUN:-0}"

CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-paper-pipeline-r}"

first_existing() {
    local path
    for path in "$@"; do
        if [ -e "${path}" ]; then
            printf '%s\n' "${path}"
            return 0
        fi
    done
    return 1
}

find_existing_output() {
    local dir="$1"
    local prefix="$2"
    shift 2

    local candidate
    for candidate in "$@"; do
        if [ -f "${dir}/${prefix}${candidate}" ]; then
            printf '%s\n' "${dir}/${prefix}${candidate}"
            return 0
        fi
    done
    return 1
}

require_path() {
    local label="$1"
    local path="$2"
    if [ ! -e "${path}" ]; then
        echo "${label} not found: ${path}" >&2
        exit 1
    fi
}

trait_file_from_source() {
    local source_file="$1"
    local trait_stem="${source_file}"
    trait_stem="${trait_stem%.summary_statistics.csv}"
    trait_stem="${trait_stem%.per_gene_estimates.tsv}"
    trait_stem="${trait_stem%.tsv.gz}"
    trait_stem="${trait_stem%.txt.gz}"
    trait_stem="${trait_stem%.csv}"
    printf '%s.per_gene_estimates.tsv\n' "${trait_stem}"
}

safe_id() {
    printf '%s' "$1" | tr -c '[:alnum:]_.-' '_'
}

ESSENTIAL_GEP="${ESSENTIAL_GEP:-$(first_existing \
    "${RUN_OUTPUT_ROOT}/perturbseq/cnmf_essential/cNMF/${ESSENTIAL_CELL}/test1/test1.gene_spectra_score.k_${K}.dt_0_4.txt" \
    "${CODE_OUTPUT_ROOT}/perturbseq/cnmf_essential/cNMF/${ESSENTIAL_CELL}/test1/test1.gene_spectra_score.k_${K}.dt_0_4.txt" || true \
)}"
ESSENTIAL_REG_DIR="${ESSENTIAL_REG_DIR:-$(first_existing \
    "${RUN_OUTPUT_ROOT}/perturbseq/cnmf_essential/cNMF_regulation/${ESSENTIAL_CELL}" \
    "${CODE_OUTPUT_ROOT}/perturbseq/cnmf_essential/cNMF_regulation/${ESSENTIAL_CELL}" || true \
)}"
GENOMEWIDE_GEP="${GENOMEWIDE_GEP:-$(first_existing \
    "${RUN_OUTPUT_ROOT}/perturbseq/cnmf_genomewide/cNMF/cNMF_all/cNMF_all.gene_spectra_score.k_${K}.dt_0_5.txt" \
    "${CODE_OUTPUT_ROOT}/perturbseq/cnmf_genomewide/cNMF/cNMF_all/cNMF_all.gene_spectra_score.k_${K}.dt_0_5.txt" || true \
)}"
GENOMEWIDE_REG_DIR="${GENOMEWIDE_REG_DIR:-$(first_existing \
    "${RUN_OUTPUT_ROOT}/perturbseq/cnmf_genomewide/cNMF_regulation/${GENOMEWIDE_CELL}" \
    "${CODE_OUTPUT_ROOT}/perturbseq/cnmf_genomewide/cNMF_regulation/${GENOMEWIDE_CELL}" || true \
)}"

require_path "MAP" "${MAP}"
require_path "LOF_DIR" "${LOF_DIR}"
require_path "GENE_MAP" "${GENE_MAP}"
require_path "SHET_PATH" "${SHET_PATH}"
require_path "R_SCRIPT" "${R_SCRIPT}"
require_path "ESSENTIAL_GEP" "${ESSENTIAL_GEP}"
require_path "ESSENTIAL_REG_DIR" "${ESSENTIAL_REG_DIR}"
require_path "GENOMEWIDE_GEP" "${GENOMEWIDE_GEP}"
require_path "GENOMEWIDE_REG_DIR" "${GENOMEWIDE_REG_DIR}"

mkdir -p \
    "${ESSENTIAL_ASSOC_DIR}" \
    "${GENOMEWIDE_ASSOC_DIR}" \
    "${JOBS_DIR}" \
    "${LOGS_DIR}" \
    "${STATUS_DIR}"

printf 'dataset\tid2\ttrait_file\tscript_path\n' > "${MANIFEST_PATH}"

generated=0
submitted=0
skipped=0
failed_submit=0

create_job() {
    local dataset="$1"
    local id2="$2"
    local trait_file="$3"
    local lof_path="$4"
    local gep_path="$5"
    local regulation_dir="$6"
    local association_dir="$7"

    local job_key script_path out_reg out_prog
    job_key="$(safe_id "${dataset}_${id2}")"
    script_path="${JOBS_DIR}/${job_key}.sbatch"
    out_reg="${association_dir}/regulators_enrichment_K${K}_${trait_file}"
    out_prog="${association_dir}/programs_enrichment_K${K}_${trait_file}"

    cat > "${script_path}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${LOGS_DIR}/${job_key}_%j.out
#SBATCH --error=${LOGS_DIR}/${job_key}_%j.err
#SBATCH --mem=${SBATCH_MEM}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${SBATCH_CPUS}
#SBATCH --partition=${SBATCH_PARTITION}
#SBATCH --time=${SBATCH_TIME}
#SBATCH --chdir=${PROJECT_ROOT}

set -euo pipefail

if [ -f "${CONDA_SH}" ]; then
    # shellcheck disable=SC1090
    source "${CONDA_SH}"
fi

export OMP_NUM_THREADS=${SBATCH_CPUS}
export MKL_NUM_THREADS=${SBATCH_CPUS}
export OPENBLAS_NUM_THREADS=${SBATCH_CPUS}
export VECLIB_MAXIMUM_THREADS=${SBATCH_CPUS}
export NUMEXPR_NUM_THREADS=${SBATCH_CPUS}

mkdir -p "${association_dir}" "${STATUS_DIR}"
printf 'time=%s\ndataset=%s\nid2=%s\njob_id=%s\n' "\$(date -Iseconds)" "${dataset}" "${id2}" "\${SLURM_JOB_ID:-local}" > "${STATUS_DIR}/${job_key}.running"

conda run --no-capture-output -n "${CONDA_ENV}" Rscript \
    "${R_SCRIPT}" \
    "${gep_path}" \
    "${GENE_MAP}" \
    "${regulation_dir}" \
    "${lof_path}" \
    "${SHET_PATH}" \
    "${K}" \
    "${out_reg}" \
    "${out_prog}" \
    "${RANDOM_ITERATIONS}" \
    "${TOP_N}"

printf 'time=%s\ndataset=%s\nid2=%s\njob_id=%s\nout_reg=%s\nout_prog=%s\n' \
    "\$(date -Iseconds)" "${dataset}" "${id2}" "\${SLURM_JOB_ID:-local}" "${out_reg}" "${out_prog}" > "${STATUS_DIR}/${job_key}.ok"
rm -f "${STATUS_DIR}/${job_key}.running" "${STATUS_DIR}/${job_key}.failed"
EOF

    chmod +x "${script_path}"
    printf '%s\t%s\t%s\t%s\n' "${dataset}" "${id2}" "${trait_file}" "${script_path}" >> "${MANIFEST_PATH}"
    generated=$((generated + 1))

    if [ "${DRY_RUN}" = "1" ]; then
        echo "[DRY RUN] sbatch ${script_path}"
        submitted=$((submitted + 1))
        return 0
    fi

    local submit_output
    if submit_output="$(sbatch --parsable "${script_path}" 2>&1)"; then
        echo "${submit_output}" > "${STATUS_DIR}/${job_key}.submitted"
        echo "[OK] ${dataset} ${id2} -> Job ${submit_output}"
        submitted=$((submitted + 1))
    else
        echo "[FAIL] ${dataset} ${id2}: ${submit_output}" >&2
        printf 'time=%s\ndataset=%s\nid2=%s\nreason=submit_failed\nmessage=%s\n' \
            "$(date -Iseconds)" "${dataset}" "${id2}" "${submit_output}" > "${STATUS_DIR}/${job_key}.failed"
        failed_submit=$((failed_submit + 1))
    fi
}

while IFS=$'\t' read -r id1 id2 path1 path2; do
    path2="${path2%$'\r'}"
    source_file="$(basename "${path2}")"
    trait_file="$(trait_file_from_source "${source_file}")"
    source_posterior_file="${source_file}.per_gene_estimates.tsv"
    legacy_trait_file="${source_file}"

    lof_path="$(first_existing \
        "${LOF_DIR}/${trait_file}" \
        "${LOF_DIR}/${source_posterior_file}" \
        "${LOF_DIR}/${legacy_trait_file}" \
    )" || {
        echo "[WARN] missing LoF posterior for ${id2}: ${trait_file}" >&2
        skipped=$((skipped + 1))
        continue
    }

    for dataset in cnmf_essential cnmf_genomewide; do
        if [ "${dataset}" = "cnmf_essential" ]; then
            assoc_dir="${ESSENTIAL_ASSOC_DIR}"
            gep_path="${ESSENTIAL_GEP}"
            reg_dir="${ESSENTIAL_REG_DIR}"
        else
            assoc_dir="${GENOMEWIDE_ASSOC_DIR}"
            gep_path="${GENOMEWIDE_GEP}"
            reg_dir="${GENOMEWIDE_REG_DIR}"
        fi

        programs_missing=0
        regulators_missing=0
        find_existing_output "${assoc_dir}" "programs_enrichment_K${K}_" "${trait_file}" "${source_posterior_file}" "${legacy_trait_file}" >/dev/null || programs_missing=1
        find_existing_output "${assoc_dir}" "regulators_enrichment_K${K}_" "${trait_file}" "${source_posterior_file}" "${legacy_trait_file}" >/dev/null || regulators_missing=1

        if [ "${programs_missing}" -eq 0 ] && [ "${regulators_missing}" -eq 0 ]; then
            skipped=$((skipped + 1))
            continue
        fi

        create_job "${dataset}" "${id2}" "${trait_file}" "${lof_path}" "${gep_path}" "${reg_dir}" "${assoc_dir}"
    done
done < <(tail -n +2 "${MAP}")

echo ""
echo "Generated sbatch scripts: ${generated}"
echo "Submitted: ${submitted}"
echo "Skipped existing/missing-input: ${skipped}"
echo "Submit failed: ${failed_submit}"
echo "Manifest: ${MANIFEST_PATH}"
echo "Logs: ${LOGS_DIR}"
echo "Status: ${STATUS_DIR}"
