#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-${SCRIPT_DIR}}"
CNMF_ALL_DIR="${CNMF_ALL_DIR:-${PROJECT_ROOT}/outputs/perturbseq/kolf/cNMF/cNMF_all}"
METADATA_PATH="${METADATA_PATH:-${PROJECT_ROOT}/data/KOLF/pipeline_ready/KOLF.cell_metadata.pipeline.csv}"
BATCH_ROOT="${BATCH_ROOT:-${PROJECT_ROOT}/scripts/kolf}"
JOBS_ROOT="${JOBS_ROOT:-${BATCH_ROOT}/jobs}"
REGULATION_DIR="${REGULATION_DIR:-${PROJECT_ROOT}/outputs/perturbseq/kolf/cNMF_regulation}"
K_VALUES="${K_VALUES:-30 60 90 120 150}"
DENSITY_LABEL="${DENSITY_LABEL:-0_5}"

MODE="check"
MAX_JOBS=0
ACTIVE_TMP=""

usage() {
    cat <<'EOF'
Usage: submit_kolf_regulation_only.sh [--check-only | --fix-only | --submit] [--max-jobs N]

  --check-only  Validate inputs without changing files or submitting jobs (default).
  --fix-only    Atomically add missing cell_barcode headers, then validate.
  --submit      Fix headers, validate, and submit only incomplete regulation jobs.
  --max-jobs N  Submit at most N jobs; use 1 for a pilot run (0 means unlimited).
EOF
}

while (($# > 0)); do
    case "$1" in
        --check-only)
            MODE="check"
            ;;
        --fix-only)
            MODE="fix"
            ;;
        --submit)
            MODE="submit"
            ;;
        --max-jobs)
            shift
            if (($# == 0)) || [[ ! "$1" =~ ^[0-9]+$ ]]; then
                echo "--max-jobs requires a non-negative integer" >&2
                exit 2
            fi
            MAX_JOBS="$1"
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
    shift
done

cleanup() {
    if [[ -n "${ACTIVE_TMP}" && -e "${ACTIVE_TMP}" ]]; then
        rm -f -- "${ACTIVE_TMP}"
    fi
}
trap cleanup EXIT INT TERM

require_file() {
    local path="$1"
    if [[ ! -s "${path}" ]]; then
        echo "Missing or empty file: ${path}" >&2
        exit 1
    fi
}

validate_metadata() {
    require_file "${METADATA_PATH}"

    if ! head -n 1 "${METADATA_PATH}" | awk -F',' '
        {
            for (i = 1; i <= NF; i++) present[$i] = 1
            required[1] = "cell_barcode"
            required[2] = "gene"
            required[3] = "gem_group"
            required[4] = "n_genes"
            required[5] = "mitopercent"
            for (i = 1; i <= 5; i++) {
                if (!(required[i] in present)) exit 1
            }
        }
    '; then
        echo "Metadata is missing a required column: ${METADATA_PATH}" >&2
        exit 1
    fi

    if ! grep -m 1 -q ',non-targeting,' "${METADATA_PATH}"; then
        echo "Metadata has no non-targeting control cells: ${METADATA_PATH}" >&2
        exit 1
    fi
}

validate_usage_file() {
    local path="$1"
    local k="$2"

    require_file "${path}"

    if ! head -n 1 "${path}" | awk -F'\t' -v k="${k}" '
        {
            if (NF != k + 1 || $1 != "cell_barcode") exit 1
            for (i = 1; i <= k; i++) {
                if ($(i + 1) != i) exit 1
            }
        }
    '; then
        echo "Invalid usages header for K=${k}: ${path}" >&2
        exit 1
    fi

    local first_fields
    local last_fields
    first_fields="$(sed -n '2p' "${path}" | awk -F'\t' '{print NF}')"
    last_fields="$(tail -n 1 "${path}" | awk -F'\t' '{print NF}')"
    if [[ "${first_fields}" != "$((k + 1))" || "${last_fields}" != "$((k + 1))" ]]; then
        echo "Usages data field count does not match K=${k}: ${path}" >&2
        exit 1
    fi
}

fix_usage_header() {
    local path="$1"
    local k="$2"
    local header
    header="$(head -n 1 "${path}")"

    if [[ "${header}" == cell_barcode$'\t'* ]]; then
        validate_usage_file "${path}" "${k}"
        echo "Header already valid: $(basename "${path}")"
        return
    fi

    if [[ "${header}" != $'\t'* ]]; then
        echo "Refusing to modify unexpected usages header: ${path}" >&2
        exit 1
    fi

    if [[ "${MODE}" == "check" ]]; then
        echo "Missing cell_barcode header: ${path}" >&2
        return 2
    fi

    local file_size
    local available_kb
    local required_kb
    file_size="$(stat -c '%s' "${path}")"
    available_kb="$(df -Pk "$(dirname "${path}")" | awk 'NR == 2 {print $4}')"
    required_kb=$((file_size / 1024 + 1048576))
    if ((available_kb < required_kb)); then
        echo "Insufficient free space to safely rewrite ${path}" >&2
        echo "Required approximately ${required_kb} KiB; available ${available_kb} KiB" >&2
        exit 1
    fi

    ACTIVE_TMP="${path}.headerfix.$$"
    {
        printf 'cell_barcode'
        cat -- "${path}"
    } > "${ACTIVE_TMP}"
    chmod --reference="${path}" "${ACTIVE_TMP}"

    local new_size
    new_size="$(stat -c '%s' "${ACTIVE_TMP}")"
    if ((new_size != file_size + 12)); then
        echo "Header rewrite size check failed for ${path}" >&2
        exit 1
    fi

    validate_usage_file "${ACTIVE_TMP}" "${k}"
    mv -f -- "${ACTIVE_TMP}" "${path}"
    ACTIVE_TMP=""
    echo "Added cell_barcode header: $(basename "${path}")"
}

regulation_output_valid() {
    local path="$1"
    [[ -s "${path}" ]] || return 1
    awk -F'\t' '
        NR == 1 {
            if ($1 != "GENE" || $2 != "lm_es" || $3 != "lm_p") {
                invalid = 1
                exit
            }
            next
        }
        NF != 3 {
            invalid = 1
            exit
        }
        {data_rows++}
        END {exit(!invalid && data_rows > 0 ? 0 : 1)}
    '
}

validate_metadata
read -r -a KS <<< "${K_VALUES}"

metadata_first="$(sed -n '2{s/,.*//;p;}' "${METADATA_PATH}")"
metadata_last="$(tail -n 1 "${METADATA_PATH}" | cut -d',' -f1)"
header_check_failed=0

for k in "${KS[@]}"; do
    if [[ ! "${k}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Invalid K value: ${k}" >&2
        exit 1
    fi

    usage_path="${CNMF_ALL_DIR}/cNMF_all.usages.k_${k}.dt_${DENSITY_LABEL}.consensus.txt"
    require_file "${usage_path}"
    if ! fix_usage_header "${usage_path}" "${k}"; then
        header_check_failed=1
        continue
    fi

    usage_first="$(sed -n '2p' "${usage_path}" | cut -f1)"
    usage_last="$(tail -n 1 "${usage_path}" | cut -f1)"
    if [[ "${usage_first}" != "${metadata_first}" || "${usage_last}" != "${metadata_last}" ]]; then
        echo "Barcode boundary mismatch between metadata and K=${k} usages" >&2
        exit 1
    fi
done

if ((header_check_failed)); then
    echo "Header validation failed. Run with --fix-only or --submit to repair the files." >&2
    exit 2
fi

if [[ ! -d "${JOBS_ROOT}" ]]; then
    echo "Missing generated jobs directory: ${JOBS_ROOT}" >&2
    exit 1
fi

declare -A COVERED=()
declare -A ACTIVE_JOBS=()
declare -A ALLOWED_K=()
declare -a ALL_SCRIPTS=()
declare -a TO_SUBMIT=()

for k in "${KS[@]}"; do
    ALLOWED_K["${k}"]=1
done

if command -v squeue >/dev/null 2>&1; then
    while IFS= read -r job_name; do
        [[ -n "${job_name}" ]] && ACTIVE_JOBS["${job_name}"]=1
    done < <(squeue -h -u "${USER}" -o '%j')

    for job_name in "${!ACTIVE_JOBS[@]}"; do
        case "${job_name}" in
            kolf_standardize|kolf_kselect_prepare|kolf_ks_fact_*|kolf_kselect_postbase|kolf_cons_*)
                echo "Upstream KOLF job is still active: ${job_name}" >&2
                exit 1
                ;;
        esac
    done
fi

mapfile -t ALL_SCRIPTS < <(
    find "${JOBS_ROOT}" -maxdepth 1 -type f -name 'regulation_K*_program_*.sh' -print | sort -V
)
if ((${#ALL_SCRIPTS[@]} == 0)); then
    echo "No regulation job scripts found under: ${JOBS_ROOT}" >&2
    exit 1
fi

for script_path in "${ALL_SCRIPTS[@]}"; do
    base="$(basename "${script_path}")"
    if [[ ! "${base}" =~ ^regulation_K([0-9]+)_program_([0-9]+)_([0-9]+)\.sh$ ]]; then
        echo "Unexpected regulation script name: ${base}" >&2
        exit 1
    fi

    k="${BASH_REMATCH[1]}"
    start="${BASH_REMATCH[2]}"
    end="${BASH_REMATCH[3]}"
    if [[ -z "${ALLOWED_K[${k}]+x}" ]]; then
        echo "Unexpected K=${k} regulation script: ${base}" >&2
        exit 1
    fi
    if ((start < 1 || end < start || end > k)); then
        echo "Invalid program range in ${base}" >&2
        exit 1
    fi

    for ((program = start; program <= end; program++)); do
        key="${k}:${program}"
        if [[ -n "${COVERED[${key}]+x}" ]]; then
            echo "Duplicate program coverage for K=${k}, program=${program}" >&2
            exit 1
        fi
        COVERED["${key}"]=1
    done
done

for k in "${KS[@]}"; do
    for ((program = 1; program <= k; program++)); do
        key="${k}:${program}"
        if [[ -z "${COVERED[${key}]+x}" ]]; then
            echo "Missing job coverage for K=${k}, program=${program}" >&2
            exit 1
        fi
    done
done

for script_path in "${ALL_SCRIPTS[@]}"; do
    base="$(basename "${script_path}")"
    [[ "${base}" =~ ^regulation_K([0-9]+)_program_([0-9]+)_([0-9]+)\.sh$ ]]
    k="${BASH_REMATCH[1]}"
    start="${BASH_REMATCH[2]}"
    end="${BASH_REMATCH[3]}"

    complete=1
    for ((program = start; program <= end; program++)); do
        output_path="${REGULATION_DIR}/K${k}_program${program}_perturb_effects.txt"
        if ! regulation_output_valid "${output_path}"; then
            complete=0
            break
        fi
    done
    if ((complete)); then
        echo "Skip completed range: K=${k}, programs=${start}-${end}"
        continue
    fi

    job_name="$(awk -F= '/^#SBATCH --job-name=/{print $2; exit}' "${script_path}")"
    if [[ -n "${job_name}" && -n "${ACTIVE_JOBS[${job_name}]+x}" ]]; then
        echo "Skip active job: ${job_name}"
        continue
    fi
    if ((MAX_JOBS == 0 || ${#TO_SUBMIT[@]} < MAX_JOBS)); then
        TO_SUBMIT+=("${script_path}")
    fi
done

echo "Validated ${#ALL_SCRIPTS[@]} regulation scripts."
echo "Incomplete ranges ready to submit: ${#TO_SUBMIT[@]}"

if [[ "${MODE}" != "submit" ]]; then
    echo "No jobs submitted (mode=${MODE})."
    exit 0
fi

if ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch is not available" >&2
    exit 1
fi

mkdir -p "${BATCH_ROOT}"
manifest="${BATCH_ROOT}/regulation_submission_$(date '+%Y%m%d_%H%M%S').tsv"
printf 'job_id\tscript\n' > "${manifest}"

for script_path in "${TO_SUBMIT[@]}"; do
    job_id="$(sbatch --parsable "${script_path}")"
    printf '%s\t%s\n' "${job_id}" "${script_path}" >> "${manifest}"
    echo "Submitted ${job_id}: $(basename "${script_path}")"
done

echo "Submission manifest: ${manifest}"
echo "Only regulation jobs were submitted; no cNMF upstream stage was rerun."
