#!/usr/bin/env bash
# 共享函数库：所有 gen_*.sh 脚本 source 此文件

set -euo pipefail

# ============================================================
# 路径默认值（与 run_figure_all.sh 一致，运行前可通过 env 覆盖）
# ============================================================
PROJECT_ROOT="${PROJECT_ROOT:-/gpfs/chencao/qinminzhang/Nature/mine/code}"

# 图和数据产出目录（运行目录下的 outputs/）
ARTIFACT_ROOT="${ARTIFACT_ROOT:-$(pwd)/outputs}"
# sbatch 脚本、提交器和日志目录
BATCH_ROOT="${BATCH_ROOT:-$(pwd)/scripts}"
LOGS_ROOT="${LOGS_ROOT:-$(pwd)/logs}"

# 输入数据路径
# cNMF regulation 数据（跑一次的，在 Nature/mine/code/outputs/ 下）
CNMF_REGULATION_DIR="${CNMF_REGULATION_DIR:-/gpfs/chencao/qinminzhang/Nature/mine/code/outputs/perturbseq/cnmf_genomewide/cNMF_regulation/K562GW}"
# cNMF spectra 矩阵（跑一次的）
SPECTRA_PATH="${SPECTRA_PATH:-/gpfs/chencao/qinminzhang/Nature/mine/code/outputs/perturbseq/cnmf_genomewide/cNMF/cNMF_all/cNMF_all.gene_spectra_score.k_60.dt_0_5.txt}"
# cNMF trait association 数据（分 trait 跑的，在 run_all/outputs/ 下）
PROGRAM_ASSOC_DIR="${PROGRAM_ASSOC_DIR:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/run_all/outputs/perturbseq/cnmf_genomewide/trait_association/K562GW/ProgramLevel}"

# posterior 在 run_all/outputs/ 下
POSTERIOR_DIR="${POSTERIOR_DIR:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/run_all/outputs/genebayes/posterior}"
# limma 矩阵 + shet 协变量（gene_level_scatter 自己算 correlation，不依赖 plot workflow）
LIMMA_PATH="${LIMMA_PATH:-/gpfs/chencao/qinminzhang/Nature/mine/code/outputs/perturbseq/gene_level/K562GW/limma_logFC_sum.txt}"
SHET_PATH="${SHET_PATH:-/gpfs/chencao/qinminzhang/Nature/mine/code/data/shet_10bins.txt}"

# 静态参考数据（硬编码绝对路径，避免 PROJECT_ROOT 被覆盖导致解析错误）
GENE_MAP="${GENE_MAP:-/gpfs/chencao/qinminzhang/Nature/mine/code/data/gencode_v41_gname_gid_ALL_sorted_onlyID}"
GENESET_DIR="${GENESET_DIR:-/gpfs/chencao/qinminzhang/Nature/mine/code/data/geneset}"
GENE_ANNOTATION="${GENE_ANNOTATION:-/gpfs/chencao/qinminzhang/Nature/mine/code/data/GWAS/genes.protein_coding.v39.gtf}"
FILE_ID_MAP="${FILE_ID_MAP:-/gpfs/chencao/qinminzhang/Nature/mine/code/configs/path.file_id_map.tsv}"
POSTERIOR_NAME_MAP="${POSTERIOR_NAME_MAP:-}"

# 图上高亮的基因集（4 个，保持图可读性）
HIGHLIGHT_GENESETS="${HIGHLIGHT_GENESETS:-HALLMARK_HEME_METABOLISM,Hematopoiesisgenes,mitotic_cell_cycle,positive_macromolecule_synthesis}"

# 基因集：默认自动扫描，也可手动指定
if [[ -z "${GENESET_LIST:-}" && -d "${GENESET_DIR}" ]]; then
    GENESET_LIST=""
    while IFS= read -r -d '' f; do
        name="$(basename "${f}" .txt)"
        [[ -z "${name}" ]] && continue
        GENESET_LIST="${GENESET_LIST}${name},"
    done < <(find "${GENESET_DIR}" -maxdepth 1 -name "*.txt" -print0 | sort -z)
    GENESET_LIST="${GENESET_LIST%,}"
fi

# Slurm / Conda 环境
CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"

# ============================================================
# 输入文件检查（每个 draw 脚本调用前先跑一遍）
# ============================================================
check_inputs() {
    local fig_name="$1"; shift
    local missing=()

    while [[ $# -gt 0 ]]; do
        local label="$1" path="$2" cat="$3" required="$4"
        shift 4
        if [[ ! -e "$path" ]]; then
            if [[ "$required" == "yes" ]]; then
                missing+=("  ${label}: ${path}  [${cat}]  ⛔ MISSING")
            else
                echo "  ${label}: ${path}  [${cat}]  ⚠ optional, not found"
            fi
        else
            echo "  ${label}: ${path}  [${cat}]  ✓"
        fi
    done

    if [[ ${#missing[@]} -gt 0 ]]; then
        echo ""
        echo "[${fig_name}] ❌ Missing required inputs:" >&2
        for m in "${missing[@]}"; do echo "$m" >&2; done
        echo "  提示: 用 env 覆盖路径，例如: ${fig_name}_LIMMA=/other/path bash cnmf.sh" >&2
        return 1
    fi
    echo ""
    return 0
}
FIGURES_R_ENV="${FIGURES_R_ENV:-paper-pipeline-plot}"
CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"

# 资源默认值
DEFAULT_MEM="${DEFAULT_MEM:-12G}"
DEFAULT_CPUS="${DEFAULT_CPUS:-2}"
DEFAULT_PARTITION="${DEFAULT_PARTITION:-cu,privority,batch01}"
DEFAULT_TIME="${DEFAULT_TIME:-12:00:00}"

# ============================================================
# 工具函数
# ============================================================
trim() { local v="$1"; v="${v#"${v%%[![:space:]]*}"}"; v="${v%"${v##*[![:space:]]}"}"; printf '%s' "$v"; }
yaml_quote() { local v="$1"; v="${v//\'/\'\'}"; printf "'%s'" "$v"; }
safe_name() { local v; v="$(trim "$1")"; v="${v//[^A-Za-z0-9._-]/_}"; printf '%s' "$v"; }

is_truthy() {
    local raw; raw="$(trim "$1")"
    case "${raw,,}" in 1|true|yes|on) return 0 ;; 0|false|no|off|'') return 1 ;; *)
        echo "Invalid boolean: ${raw}" >&2; exit 1 ;; esac
}

resolve_path() {
    local raw="$1"
    [[ "$raw" = /* ]] && { printf '%s' "$raw"; return; }
    (cd "${PROJECT_ROOT}" && realpath -m "$raw")
}

derive_posterior_basename() {
    local path="$1"; local stem; stem="$(basename "$path")"
    for s in .summary_statistics.csv .csv.gz .tsv.gz .txt.gz .csv .tsv; do
        [[ "$stem" == *"$s" ]] && { stem="${stem%$s}"; break; }
    done
    printf '%s.per_gene_estimates.tsv' "$stem"
}

resolve_trait_file() {
    local lof_id="$1"
    # 先查 name_map
    if [[ -n "${POSTERIOR_NAME_MAP}" && -f "${POSTERIOR_NAME_MAP}" ]]; then
        local fallback; fallback="$(awk -F'\t' -v id="$lof_id" '$1==id{print $2; exit}' "${POSTERIOR_NAME_MAP}")"
        [[ -n "$fallback" ]] && { printf '%s' "$fallback"; return; }
    fi
    # 从 file_id_map 推导
    local path2; path2="$(awk -F'\t' -v id="$lof_id" 'NR>1 && $2==id {print $4; exit}' "${FILE_ID_MAP}")"
    [[ -z "$path2" ]] && { echo "Cannot resolve trait file for $lof_id" >&2; exit 1; }
    derive_posterior_basename "$path2"
}

batch_id() {
    local prefix="$1"; shift; local ids=("$@")
    printf '%s_n%s__%s__%s' "$prefix" "${#ids[@]}" "${ids[0]}" "${ids[-1]}"
}

# ============================================================
# YAML 配置生成
# ============================================================
write_config_header() {
    local cfg="$1"
    cat > "$cfg" <<YEOF
project_root: $(yaml_quote "$PROJECT_ROOT")
artifact_root: $(yaml_quote "$ARTIFACT_ROOT")

executables:
  plot_rscript:
    - conda
    - run
    - --no-capture-output
    - -n
    - ${FIGURES_R_ENV}
    - Rscript

workflows:
  figures:
    shared_inputs:
      file_id_map: $(yaml_quote "$FILE_ID_MAP")
      gene_map: $(yaml_quote "$GENE_MAP")
      geneset_dir: $(yaml_quote "$GENESET_DIR")
      gene_annotation: $(yaml_quote "$GENE_ANNOTATION")
      posterior_dir: $(yaml_quote "$POSTERIOR_DIR")
YEOF
    [[ -n "$POSTERIOR_NAME_MAP" ]] && cat >> "$cfg" <<YEOF
      posterior_name_map: $(yaml_quote "$POSTERIOR_NAME_MAP")
YEOF
}

# 基因集 YAML 列表（预计算一次）
_build_gs_yaml() {
    local result=""; IFS=',' read -r -a arr <<< "${1}"
    for gs in "${arr[@]}"; do
        [[ -z "$gs" ]] && continue
        gs="${gs//\'/\'\'}"
        result+="          - '${gs}'"$'\n'
    done
    eval "$2=\"\$result\""
}

# ============================================================
# sbatch 脚本生成
# ============================================================
write_slurm_script() {
    local spath="$1" jname="$2" mem="$3" cpu="$4" part="$5" time="$6" wf="$7" cfg="$8"
    local sdir ldir
    sdir="$(cd "$(dirname "$spath")" && pwd)"
    ldir="${LOG_DIR}/${jname%%_*}"
    mkdir -p "$ldir"
    cat > "$spath" <<YEOF
#!/usr/bin/env bash
#SBATCH --job-name=${jname}
#SBATCH --error=${ldir}/${jname}.error
#SBATCH --output=${ldir}/${jname}.out
#SBATCH --mem=${mem}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpu}
#SBATCH --partition=${part}
#SBATCH --time=${time}
#SBATCH --chdir=${sdir}

set -euo pipefail
source "${CONDA_SH}"
conda activate "${CONTROL_ENV}"
cd "${sdir}"
paper-pipeline run --config "${cfg}" ${wf}
YEOF
    chmod +x "$spath"
}

write_submit() {
    local spath="$1"; shift
    { echo '#!/usr/bin/env bash'; echo 'set -euo pipefail';
      for line in "$@"; do echo "$line"; done; } > "$spath"
    chmod +x "$spath"
}
