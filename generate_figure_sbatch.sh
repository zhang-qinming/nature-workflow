#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# 路径与通用设置
# ============================================================
DEFAULT_PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BATCH_ROOT="${BATCH_ROOT:-${DEFAULT_PROJECT_ROOT}/scripts}"
CONFIG_NAME="${CONFIG_NAME:-pipeline.figures.yaml}"
LOGS_ROOT="${LOGS_ROOT:-${DEFAULT_PROJECT_ROOT}/logs/figures}"
# 并行度：同时生成多少个 trait 的脚本（xargs -P）
FIGURE_PARALLEL="${FIGURE_PARALLEL:-20}"

CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"
FIGURES_R_ENV="${FIGURES_R_ENV:-paper-pipeline-plot}"
CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"

# ============================================================
# 全局默认资源（作为各图的兜底值）
# 可用分区: cu,privority,batch01 (56C/256G), fat (56C/1T)
# ============================================================
FIGURE_DEFAULT_MEM="${FIGURE_DEFAULT_MEM:-12G}"
FIGURE_DEFAULT_CPUS="${FIGURE_DEFAULT_CPUS:-2}"
FIGURE_DEFAULT_PARTITION="${FIGURE_DEFAULT_PARTITION:-cu,privority,batch01,fat}"
FIGURE_DEFAULT_TIME="${FIGURE_DEFAULT_TIME:-12:00:00}"

# ============================================================
# 每种图的资源量（mem / cpus / partition / time 独立控制）
# 未设置时回退到 FIGURE_DEFAULT_* 兜底值
# ============================================================

# --- cNMF 图：program-regulator 散点 + corregulation 散点 ---
FIGURE_CNMF_MEM="${FIGURE_CNMF_MEM:-8G}"
FIGURE_CNMF_CPUS="${FIGURE_CNMF_CPUS:-${FIGURE_DEFAULT_CPUS}}"
FIGURE_CNMF_PARTITION="${FIGURE_CNMF_PARTITION:-${FIGURE_DEFAULT_PARTITION}}"
FIGURE_CNMF_TIME="${FIGURE_CNMF_TIME:-${FIGURE_DEFAULT_TIME}}"

# --- LoF burden 火山图 ---
FIGURE_BURDEN_VOLCANO_MEM="${FIGURE_BURDEN_VOLCANO_MEM:-8G}"
FIGURE_BURDEN_VOLCANO_CPUS="${FIGURE_BURDEN_VOLCANO_CPUS:-${FIGURE_DEFAULT_CPUS}}"
FIGURE_BURDEN_VOLCANO_PARTITION="${FIGURE_BURDEN_VOLCANO_PARTITION:-${FIGURE_DEFAULT_PARTITION}}"
FIGURE_BURDEN_VOLCANO_TIME="${FIGURE_BURDEN_VOLCANO_TIME:-${FIGURE_DEFAULT_TIME}}"

# --- 基因级散点图 (posterior vs perturb-seq) ---
FIGURE_GENE_LEVEL_SCATTER_MEM="${FIGURE_GENE_LEVEL_SCATTER_MEM:-8G}"
FIGURE_GENE_LEVEL_SCATTER_CPUS="${FIGURE_GENE_LEVEL_SCATTER_CPUS:-${FIGURE_DEFAULT_CPUS}}"
FIGURE_GENE_LEVEL_SCATTER_PARTITION="${FIGURE_GENE_LEVEL_SCATTER_PARTITION:-${FIGURE_DEFAULT_PARTITION}}"
FIGURE_GENE_LEVEL_SCATTER_TIME="${FIGURE_GENE_LEVEL_SCATTER_TIME:-${FIGURE_DEFAULT_TIME}}"

# --- Program 排名条形图 ---
FIGURE_PROGRAM_RANKINGS_MEM="${FIGURE_PROGRAM_RANKINGS_MEM:-4G}"
FIGURE_PROGRAM_RANKINGS_CPUS="${FIGURE_PROGRAM_RANKINGS_CPUS:-1}"
FIGURE_PROGRAM_RANKINGS_PARTITION="${FIGURE_PROGRAM_RANKINGS_PARTITION:-${FIGURE_DEFAULT_PARTITION}}"
FIGURE_PROGRAM_RANKINGS_TIME="${FIGURE_PROGRAM_RANKINGS_TIME:-${FIGURE_DEFAULT_TIME}}"

# --- GWAS Manhattan 图（★重, 全基因组变异位点 + 54 基因集注释 + 光栅化） ---
FIGURE_GWAS_MANHATTAN_MEM="${FIGURE_GWAS_MANHATTAN_MEM:-64G}"
FIGURE_GWAS_MANHATTAN_CPUS="${FIGURE_GWAS_MANHATTAN_CPUS:-4}"
FIGURE_GWAS_MANHATTAN_PARTITION="${FIGURE_GWAS_MANHATTAN_PARTITION:-${FIGURE_DEFAULT_PARTITION}}"
FIGURE_GWAS_MANHATTAN_TIME="${FIGURE_GWAS_MANHATTAN_TIME:-24:00:00}"

# --- Program × Trait 热图（★重, N trait 合并 + hclust 聚类） ---
FIGURE_PROGRAM_HEATMAP_MEM="${FIGURE_PROGRAM_HEATMAP_MEM:-32G}"
FIGURE_PROGRAM_HEATMAP_CPUS="${FIGURE_PROGRAM_HEATMAP_CPUS:-4}"
FIGURE_PROGRAM_HEATMAP_PARTITION="${FIGURE_PROGRAM_HEATMAP_PARTITION:-${FIGURE_DEFAULT_PARTITION}}"
FIGURE_PROGRAM_HEATMAP_TIME="${FIGURE_PROGRAM_HEATMAP_TIME:-${FIGURE_DEFAULT_TIME}}"

# --- 基因级 QQ 图（★重, N trait × ~20000 基因合并） ---
FIGURE_GENE_LEVEL_QQ_MEM="${FIGURE_GENE_LEVEL_QQ_MEM:-32G}"
FIGURE_GENE_LEVEL_QQ_CPUS="${FIGURE_GENE_LEVEL_QQ_CPUS:-4}"
FIGURE_GENE_LEVEL_QQ_PARTITION="${FIGURE_GENE_LEVEL_QQ_PARTITION:-${FIGURE_DEFAULT_PARTITION}}"
FIGURE_GENE_LEVEL_QQ_TIME="${FIGURE_GENE_LEVEL_QQ_TIME:-${FIGURE_DEFAULT_TIME}}"

# --- Cross-Trait 相关性热图（★重, N trait posterior merge + N×N 相关矩阵） ---
FIGURE_CROSS_TRAIT_HEATMAP_MEM="${FIGURE_CROSS_TRAIT_HEATMAP_MEM:-32G}"
FIGURE_CROSS_TRAIT_HEATMAP_CPUS="${FIGURE_CROSS_TRAIT_HEATMAP_CPUS:-4}"
FIGURE_CROSS_TRAIT_HEATMAP_PARTITION="${FIGURE_CROSS_TRAIT_HEATMAP_PARTITION:-${FIGURE_DEFAULT_PARTITION}}"
FIGURE_CROSS_TRAIT_HEATMAP_TIME="${FIGURE_CROSS_TRAIT_HEATMAP_TIME:-${FIGURE_DEFAULT_TIME}}"

# --- cNMF program top 基因条形图 ---
FIGURE_CNMF_PROGRAM_TOP_GENES_MEM="${FIGURE_CNMF_PROGRAM_TOP_GENES_MEM:-8G}"
FIGURE_CNMF_PROGRAM_TOP_GENES_CPUS="${FIGURE_CNMF_PROGRAM_TOP_GENES_CPUS:-${FIGURE_DEFAULT_CPUS}}"
FIGURE_CNMF_PROGRAM_TOP_GENES_PARTITION="${FIGURE_CNMF_PROGRAM_TOP_GENES_PARTITION:-${FIGURE_DEFAULT_PARTITION}}"
FIGURE_CNMF_PROGRAM_TOP_GENES_TIME="${FIGURE_CNMF_PROGRAM_TOP_GENES_TIME:-${FIGURE_DEFAULT_TIME}}"

# --- cNMF program 基因集富集气泡图 ---
FIGURE_CNMF_PROGRAM_ENRICHMENT_MEM="${FIGURE_CNMF_PROGRAM_ENRICHMENT_MEM:-8G}"
FIGURE_CNMF_PROGRAM_ENRICHMENT_CPUS="${FIGURE_CNMF_PROGRAM_ENRICHMENT_CPUS:-${FIGURE_DEFAULT_CPUS}}"
FIGURE_CNMF_PROGRAM_ENRICHMENT_PARTITION="${FIGURE_CNMF_PROGRAM_ENRICHMENT_PARTITION:-${FIGURE_DEFAULT_PARTITION}}"
FIGURE_CNMF_PROGRAM_ENRICHMENT_TIME="${FIGURE_CNMF_PROGRAM_ENRICHMENT_TIME:-${FIGURE_DEFAULT_TIME}}"

# --- Cross-Trait 散点图（成对） ---
FIGURE_CROSS_TRAIT_MEM="${FIGURE_CROSS_TRAIT_MEM:-8G}"
FIGURE_CROSS_TRAIT_CPUS="${FIGURE_CROSS_TRAIT_CPUS:-2}"
FIGURE_CROSS_TRAIT_PARTITION="${FIGURE_CROSS_TRAIT_PARTITION:-${FIGURE_DEFAULT_PARTITION}}"
FIGURE_CROSS_TRAIT_TIME="${FIGURE_CROSS_TRAIT_TIME:-${FIGURE_DEFAULT_TIME}}"

# --- GWAS Locus Zoom ---
FIGURE_GWAS_LOCUS_ZOOM_MEM="${FIGURE_GWAS_LOCUS_ZOOM_MEM:-8G}"
FIGURE_GWAS_LOCUS_ZOOM_CPUS="${FIGURE_GWAS_LOCUS_ZOOM_CPUS:-2}"
FIGURE_GWAS_LOCUS_ZOOM_PARTITION="${FIGURE_GWAS_LOCUS_ZOOM_PARTITION:-${FIGURE_DEFAULT_PARTITION}}"
FIGURE_GWAS_LOCUS_ZOOM_TIME="${FIGURE_GWAS_LOCUS_ZOOM_TIME:-${FIGURE_DEFAULT_TIME}}"

# ============================================================

# 数据源配置
# ============================================================
# 基因-性状-文件映射表（4 列 TSV：id1, id2, path1, path2）
# 不设时自动查找 ${PROJECT_ROOT}/configs/path.file_id_map.tsv
FIGURE_FILE_ID_MAP="${FIGURE_FILE_ID_MAP:-}"
# 后验文件名映射表（可选；2 列 TSV：lof_id, posterior_filename）
FIGURE_POSTERIOR_NAME_MAP="${FIGURE_POSTERIOR_NAME_MAP:-}"
# 后验文件是否必须已存在（0=用推导名即可, 1=不存在时报错）
FIGURE_REQUIRE_EXISTING_POSTERIOR="${FIGURE_REQUIRE_EXISTING_POSTERIOR:-0}"
# 限定画图的 LoF / GWAS ID 子集（逗号分隔；不设时从 file_id_map 取全部）
FIGURE_LOF_IDS="${FIGURE_LOF_IDS:-}"
FIGURE_GWAS_IDS="${FIGURE_GWAS_IDS:-}"
# 热图 / QQ 图的 LoF ID 子集（逗号分隔；不设时用全部 lof_ids）
FIGURE_HEATMAP_LOF_IDS="${FIGURE_HEATMAP_LOF_IDS:-}"
FIGURE_QQ_LOF_IDS="${FIGURE_QQ_LOF_IDS:-}"
# Cross-Trait 成对图：id_x:id_y,id_x2:id_y2,...
FIGURE_CROSS_TRAIT_PAIRS="${FIGURE_CROSS_TRAIT_PAIRS:-}"
# Locus Zoom 的坐标配置文件（8 列 TSV）
FIGURE_LOCI_TSV="${FIGURE_LOCI_TSV:-}"
# 每种图的开关（1=启用, 0=禁用）
# ============================================================
FIGURE_ENABLE_CNMF="${FIGURE_ENABLE_CNMF:-1}"
FIGURE_ENABLE_PROGRAM_RANKINGS="${FIGURE_ENABLE_PROGRAM_RANKINGS:-1}"
FIGURE_ENABLE_PROGRAM_HEATMAP="${FIGURE_ENABLE_PROGRAM_HEATMAP:-1}"
FIGURE_ENABLE_BURDEN_VOLCANO="${FIGURE_ENABLE_BURDEN_VOLCANO:-1}"
FIGURE_ENABLE_GWAS_MANHATTAN="${FIGURE_ENABLE_GWAS_MANHATTAN:-1}"
FIGURE_ENABLE_GENE_LEVEL_SCATTER="${FIGURE_ENABLE_GENE_LEVEL_SCATTER:-1}"
FIGURE_ENABLE_GENE_LEVEL_QQ="${FIGURE_ENABLE_GENE_LEVEL_QQ:-1}"
FIGURE_ENABLE_CROSS_TRAIT="${FIGURE_ENABLE_CROSS_TRAIT:-0}"
FIGURE_ENABLE_CROSS_TRAIT_HEATMAP="${FIGURE_ENABLE_CROSS_TRAIT_HEATMAP:-1}"
FIGURE_ENABLE_GWAS_LOCUS_ZOOM="${FIGURE_ENABLE_GWAS_LOCUS_ZOOM:-auto}"
FIGURE_ENABLE_CNMF_PROGRAM_TOP_GENES="${FIGURE_ENABLE_CNMF_PROGRAM_TOP_GENES:-1}"
FIGURE_ENABLE_CNMF_PROGRAM_ENRICHMENT="${FIGURE_ENABLE_CNMF_PROGRAM_ENRICHMENT:-1}"

# ============================================================
# 图参数
# ============================================================
FIGURE_GENESET_LIST="${FIGURE_GENESET_LIST:-}"
FIGURE_CNMF_K="${FIGURE_CNMF_K:-60}"
FIGURE_PROGRAM_TOP_N="${FIGURE_PROGRAM_TOP_N:-12}"
FIGURE_GENE_LEVEL_TOP_N="${FIGURE_GENE_LEVEL_TOP_N:-8}"
FIGURE_GENE_LEVEL_Y_LIMIT="${FIGURE_GENE_LEVEL_Y_LIMIT:-8}"
FIGURE_QQ_Y_LIMIT="${FIGURE_QQ_Y_LIMIT:-10}"
FIGURE_GWAS_FLANK_BP="${FIGURE_GWAS_FLANK_BP:-50000}"
FIGURE_GWAS_LABEL_P_THRESHOLD="${FIGURE_GWAS_LABEL_P_THRESHOLD:-1e-30}"
FIGURE_GWAS_GENOMEWIDE_THRESHOLD="${FIGURE_GWAS_GENOMEWIDE_THRESHOLD:-5e-8}"
FIGURE_LOCUS_FLANK_BP="${FIGURE_LOCUS_FLANK_BP:-250000}"
FIGURE_LOCUS_LABEL_TOP_N="${FIGURE_LOCUS_LABEL_TOP_N:-6}"
FIGURE_CROSS_TRAIT_TOP_N="${FIGURE_CROSS_TRAIT_TOP_N:-12}"
FIGURE_CROSS_TRAIT_METHOD="${FIGURE_CROSS_TRAIT_METHOD:-pearson}"

# ============================================================

# ============================================================
# 路径初始化
# ============================================================
mkdir -p "${BATCH_ROOT}" "${LOGS_ROOT}"
BATCH_ROOT="$(cd "${BATCH_ROOT}" && pwd)"
LOGS_ROOT="$(cd "${LOGS_ROOT}" && pwd)"

PROJECT_ROOT="${PROJECT_ROOT:-}"
if [[ -z "${PROJECT_ROOT}" ]]; then
    PROJECT_ROOT="$(cd "${BATCH_ROOT}/../.." && pwd)"
fi
if [[ ! -d "${PROJECT_ROOT}" ]]; then
    echo "PROJECT_ROOT does not exist: ${PROJECT_ROOT}" >&2
    exit 1
fi
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${PROJECT_ROOT}/outputs}"
if [[ "${ARTIFACT_ROOT}" != /* ]]; then
    ARTIFACT_ROOT="$(cd "${PROJECT_ROOT}" && realpath -m "${ARTIFACT_ROOT}")"
fi

FIGURE_POSTERIOR_DIR="${FIGURE_POSTERIOR_DIR:-${ARTIFACT_ROOT}/genebayes/posterior}"
# cNMF regulation（跑一次）在 Nature/mine/code/outputs，association（分 trait）在 run_all/outputs
FIGURE_CNMF_REGULATION_DIR="${FIGURE_CNMF_REGULATION_DIR:-/gpfs/chencao/qinminzhang/Nature/mine/code/outputs/perturbseq/cnmf_genomewide/cNMF_regulation/K562GW}"
FIGURE_PROGRAM_ASSOCIATION_DIR="${FIGURE_PROGRAM_ASSOCIATION_DIR:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/run_all/outputs/perturbseq/cnmf_genomewide/trait_association/K562GW/ProgramLevel}"
FIGURE_LIMMA_PATH="${FIGURE_LIMMA_PATH:-/gpfs/chencao/qinminzhang/Nature/mine/code/outputs/perturbseq/gene_level/K562GW/limma_logFC_sum.txt}"
FIGURE_SHET_PATH="${FIGURE_SHET_PATH:-/gpfs/chencao/qinminzhang/Nature/mine/code/data/shet_10bins.txt}"
FIGURE_SPECTRA_PATH="${FIGURE_SPECTRA_PATH:-/gpfs/chencao/qinminzhang/Nature/mine/code/outputs/perturbseq/cnmf_genomewide/cNMF/cNMF_all/cNMF_all.gene_spectra_score.k_60.dt_0_5.txt}"
FIGURE_GENE_MAP="${FIGURE_GENE_MAP:-/gpfs/chencao/qinminzhang/Nature/mine/code/data/gencode_v41_gname_gid_ALL_sorted_onlyID}"
FIGURE_GENESET_DIR="${FIGURE_GENESET_DIR:-/gpfs/chencao/qinminzhang/Nature/mine/code/data/geneset}"
FIGURE_GENE_ANNOTATION="${FIGURE_GENE_ANNOTATION:-/gpfs/chencao/qinminzhang/Nature/mine/code/data/GWAS/genes.protein_coding.v39.gtf}"

# 输出子目录（worker 函数引用，必须在它们定义之前）
LOF_ROOT="${BATCH_ROOT}/lof"
GWAS_ROOT="${BATCH_ROOT}/gwas"
GROUP_ROOT="${BATCH_ROOT}/groups"
PAIR_ROOT="${BATCH_ROOT}/pairs"
LOCUS_ROOT="${BATCH_ROOT}/loci"

# ============================================================
# 基因集列表自动扫描
# ============================================================
# 手动指定示例：export FIGURE_GENESET_LIST="HALLMARK_HEME_METABOLISM,HALLMARK_HYPOXIA,..."
if [[ -z "${FIGURE_GENESET_LIST}" ]]; then
    if [[ -d "${FIGURE_GENESET_DIR}" ]]; then
        FIGURE_GENESET_LIST=""
        while IFS= read -r -d '' f; do
            name="$(basename "${f}" .txt)"
            [[ -z "${name}" ]] && continue
            FIGURE_GENESET_LIST="${FIGURE_GENESET_LIST}${name},"
        done < <(find "${FIGURE_GENESET_DIR}" -maxdepth 1 -name "*.txt" -print0 | sort -z)
        FIGURE_GENESET_LIST="${FIGURE_GENESET_LIST%,}"
    fi
    if [[ -z "${FIGURE_GENESET_LIST}" ]]; then
        FIGURE_GENESET_LIST="Autophagosome_genes,HALLMARK_ADIPOGENESIS,HALLMARK_ALLOGRAFT_REJECTION,HALLMARK_ANDROGEN_RESPONSE,HALLMARK_ANGIOGENESIS,HALLMARK_APICAL_JUNCTION,HALLMARK_APICAL_SURFACE,HALLMARK_APOPTOSIS,HALLMARK_BILE_ACID_METABOLISM,HALLMARK_CHOLESTEROL_HOMEOSTASIS,HALLMARK_COAGULATION,HALLMARK_COMPLEMENT,HALLMARK_DNA_REPAIR,HALLMARK_E2F_TARGETS,HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,HALLMARK_ESTROGEN_RESPONSE_EARLY,HALLMARK_ESTROGEN_RESPONSE_LATE,HALLMARK_FATTY_ACID_METABOLISM,HALLMARK_G2M_CHECKPOINT,HALLMARK_GLYCOLYSIS,HALLMARK_HEDGEHOG_SIGNALING,HALLMARK_HEME_METABOLISM,HALLMARK_HYPOXIA,HALLMARK_IL2_STAT5_SIGNALING,HALLMARK_IL6_JAK_STAT3_SIGNALING,HALLMARK_INFLAMMATORY_RESPONSE,HALLMARK_INTERFERON_ALPHA_RESPONSE,HALLMARK_INTERFERON_GAMMA_RESPONSE,HALLMARK_KRAS_SIGNALING_DN,HALLMARK_KRAS_SIGNALING_UP,HALLMARK_MITOTIC_SPINDLE,HALLMARK_MTORC1_SIGNALING,HALLMARK_MYC_TARGETS_V1,HALLMARK_MYC_TARGETS_V2,HALLMARK_MYOGENESIS,HALLMARK_NOTCH_SIGNALING,HALLMARK_OXIDATIVE_PHOSPHORYLATION,HALLMARK_P53_PATHWAY,HALLMARK_PANCREAS_BETA_CELLS,HALLMARK_PEROXISOME,HALLMARK_PI3K_AKT_MTOR_SIGNALING,HALLMARK_PROTEIN_SECRETION,HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,HALLMARK_SPERMATOGENESIS,HALLMARK_TGF_BETA_SIGNALING,HALLMARK_TNFA_SIGNALING_VIA_NFKB,HALLMARK_UNFOLDED_PROTEIN_RESPONSE,HALLMARK_UV_RESPONSE_DN,HALLMARK_UV_RESPONSE_UP,HALLMARK_WNT_BETA_CATENIN_SIGNALING,HALLMARK_XENOBIOTIC_METABOLISM,Hematopoiesisgenes,mitotic_cell_cycle,positive_macromolecule_synthesis"
    fi
fi

# ============================================================
# 工具函数
# ============================================================
trim() {
    local value="$1"
    value="${value#"${value%%[![:space:]]*}"}"
    value="${value%"${value##*[![:space:]]}"}"
    printf '%s' "$value"
}

yaml_quote() {
    local value="$1"
    value="${value//\'/\'\'}"
    printf "'%s'" "$value"
}

resolve_path() {
    local raw_path="$1"
    if [[ "$raw_path" = /* ]]; then
        printf '%s' "$raw_path"
    else
        (
            cd "${PROJECT_ROOT}"
            realpath -m "$raw_path"
        )
    fi
}

is_truthy() {
    local raw
    raw="$(trim "$1")"
    case "${raw,,}" in
        1|true|yes|on) return 0 ;;
        0|false|no|off|'') return 1 ;;
        *)
            echo "Boolean parameter must be one of: 0/1, true/false, yes/no, on/off. Got: ${raw}" >&2
            exit 1
            ;;
    esac
}

safe_name() {
    local raw="$1"
    raw="$(trim "$raw")"
    raw="${raw//[^A-Za-z0-9._-]/_}"
    printf '%s' "$raw"
}

# ============================================================
# file_id_map 校验
# ============================================================
if [[ -z "${FIGURE_FILE_ID_MAP}" ]]; then
    if [[ -f "${PROJECT_ROOT}/configs/path.file_id_map.tsv" ]]; then
        FIGURE_FILE_ID_MAP="${PROJECT_ROOT}/configs/path.file_id_map.tsv"
    else
        echo "FIGURE_FILE_ID_MAP is required." >&2
        echo "  Set it: export FIGURE_FILE_ID_MAP=/path/to/mapping.tsv" >&2
        echo "  Or create: ${PROJECT_ROOT}/configs/path.file_id_map.tsv" >&2
        echo "  Format: TSV with columns id1, id2, path1, path2" >&2
        exit 1
    fi
fi
FIGURE_FILE_ID_MAP="$(resolve_path "${FIGURE_FILE_ID_MAP}")"
if [[ ! -f "${FIGURE_FILE_ID_MAP}" ]]; then
    echo "FIGURE_FILE_ID_MAP does not exist: ${FIGURE_FILE_ID_MAP}" >&2
    exit 1
fi

IFS=$'\t' read -r _h1 _h2 _h3 _h4 _rest < "${FIGURE_FILE_ID_MAP}"
if [[ "${_h1}" != "id1" || "${_h2}" != "id2" || "${_h3}" != "path1" || "${_h4}" != "path2" ]]; then
    echo "FIGURE_FILE_ID_MAP must be a TSV with header columns: id1, id2, path1, path2" >&2
    echo "The old 2-column configs/file_id_map.tsv is not valid here." >&2
    exit 1
fi
FIGURE_FILE_ID_MAP_DIR="$(cd "$(dirname "${FIGURE_FILE_ID_MAP}")" && pwd)"

if [[ -n "${FIGURE_POSTERIOR_NAME_MAP}" ]]; then
    FIGURE_POSTERIOR_NAME_MAP="$(resolve_path "${FIGURE_POSTERIOR_NAME_MAP}")"
    if [[ ! -f "${FIGURE_POSTERIOR_NAME_MAP}" ]]; then
        echo "FIGURE_POSTERIOR_NAME_MAP does not exist: ${FIGURE_POSTERIOR_NAME_MAP}" >&2
        exit 1
    fi
fi

resolve_from_map_dir() {
    local raw_path="$1"
    if [[ "$raw_path" = /* ]]; then
        printf '%s' "$raw_path"
    else
        (
            cd "${FIGURE_FILE_ID_MAP_DIR}"
            realpath -m "$raw_path"
        )
    fi
}

map_path1_by_id() {
    local source_id="$1"
    local raw
    raw="$(awk -F'\t' -v id="${source_id}" 'NR > 1 && $1 == id {print $3; exit}' "${FIGURE_FILE_ID_MAP}")"
    [[ -n "${raw}" ]] || return 1
    resolve_from_map_dir "${raw}"
}

map_path2_by_id() {
    local source_id="$1"
    local raw
    raw="$(awk -F'\t' -v id="${source_id}" 'NR > 1 && $2 == id {print $4; exit}' "${FIGURE_FILE_ID_MAP}")"
    [[ -n "${raw}" ]] || return 1
    resolve_from_map_dir "${raw}"
}

collect_map_ids() {
    local column="$1"
    awk -F'\t' -v col="${column}" 'NR > 1 && $col != "" {print $col}' "${FIGURE_FILE_ID_MAP}" | sort -u
}

array_from_csv() {
    local raw="$1"
    local out_var="$2"
    local item
    eval "${out_var}=()"
    IFS=',' read -r -a _items <<< "${raw}"
    for item in "${_items[@]}"; do
        item="$(trim "$item")"
        [[ -z "${item}" ]] && continue
        eval "${out_var}+=(\"${item}\")"
    done
}

derive_posterior_basename_from_source() {
    local source_path="$1"
    local stem
    stem="$(basename "${source_path}")"
    for suffix in \
        ".summary_statistics.csv" \
        ".csv.gz" \
        ".tsv.gz" \
        ".txt.gz" \
        ".csv" \
        ".tsv"; do
        if [[ "${stem}" == *"${suffix}" ]]; then
            stem="${stem%${suffix}}"
            break
        fi
    done
    printf '%s.per_gene_estimates.tsv' "${stem}"
}

lookup_posterior_name_map() {
    local lof_id="$1"
    [[ -n "${FIGURE_POSTERIOR_NAME_MAP}" ]] || return 1
    awk -F'\t' -v id="${lof_id}" '
        NR == 1 && ($1 == "id" || $1 == "id2" || $1 == "lof_id") {next}
        $1 == id {print $2; exit}
    ' "${FIGURE_POSTERIOR_NAME_MAP}"
}

resolve_trait_file() {
    local lof_id="$1"
    local mapped_path derived match fallback

    fallback="$(lookup_posterior_name_map "${lof_id}" || true)"
    if [[ -n "${fallback}" ]]; then
        if [[ -f "${FIGURE_POSTERIOR_DIR}/${fallback}" ]]; then
            printf '%s' "${fallback}"
            return 0
        fi
        if is_truthy "${FIGURE_REQUIRE_EXISTING_POSTERIOR}"; then
            echo "Posterior name map resolved ${lof_id} -> ${fallback}, but file does not exist under ${FIGURE_POSTERIOR_DIR}" >&2
            exit 1
        fi
        printf '%s' "${fallback}"
        return 0
    fi

    mapped_path="$(map_path2_by_id "${lof_id}" || true)"
    if [[ -z "${mapped_path}" ]]; then
        echo "Unable to resolve path2 for LoF ID ${lof_id}" >&2
        exit 1
    fi
    derived="$(derive_posterior_basename_from_source "${mapped_path}")"
    if [[ -f "${FIGURE_POSTERIOR_DIR}/${derived}" ]]; then
        printf '%s' "${derived}"
        return 0
    fi

    shopt -s nullglob
    local matches=("${FIGURE_POSTERIOR_DIR}"/*"${lof_id}"*.per_gene_estimates.tsv)
    shopt -u nullglob
    if (( ${#matches[@]} >= 1 )); then
        basename "${matches[0]}"
        return 0
    fi

    if is_truthy "${FIGURE_REQUIRE_EXISTING_POSTERIOR}"; then
        echo "Unable to resolve posterior file for LoF ID ${lof_id}. Expected ${FIGURE_POSTERIOR_DIR}/${derived}" >&2
        exit 1
    fi
    printf '%s' "${derived}"
}

batch_id_from_ids() {
    local prefix="$1"
    shift
    local ids=("$@")
    local count="${#ids[@]}"
    local first="${ids[0]}"
    local last="${ids[$((count - 1))]}"
    printf '%s_n%s__%s__%s' "${prefix}" "${count}" "${first}" "${last}"
}

# ============================================================
# 输出生成函数
# ============================================================
write_slurm_script() {
    local script_path="$1"
    local job_name="$2"
    local mem="$3"
    local cpus="$4"
    local partition="$5"
    local time_limit="$6"
    local workflow_name="$7"
    local config_path="$8"
    local log_subdir="${9:-}"
    local script_dir log_dir

    script_dir="$(cd "$(dirname "${script_path}")" && pwd)"
    log_dir="${LOGS_ROOT}/${log_subdir}"
    mkdir -p "${log_dir}"

    cat > "${script_path}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${job_name}
#SBATCH --error=${log_dir}/${job_name}.error
#SBATCH --output=${log_dir}/${job_name}.out
#SBATCH --mem=${mem}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --partition=${partition}
#SBATCH --time=${time_limit}
#SBATCH --chdir=${script_dir}

set -euo pipefail

source "${CONDA_SH}"
conda activate "${CONTROL_ENV}"

cd "${script_dir}"
paper-pipeline run --config "${config_path}" ${workflow_name}
EOF
    chmod +x "${script_path}"
}

write_submit_helper_from_arrays() {
    local script_path="$1"
    shift
    local array_name count idx line

    {
        echo '#!/usr/bin/env bash'
        echo 'set -euo pipefail'
        for array_name in "$@"; do
            [[ "${array_name}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || {
                echo "Invalid submit array variable name: ${array_name}" >&2
                exit 1
            }
            eval "count=\${#${array_name}[@]}"
            for ((idx = 0; idx < count; idx++)); do
                eval "line=\${${array_name}[idx]}"
                echo "${line}"
            done
        done
    } > "${script_path}"
    chmod +x "${script_path}"
}

write_config_header() {
    local config_path="$1"
    cat > "${config_path}" <<EOF
project_root: $(yaml_quote "${PROJECT_ROOT}")
artifact_root: $(yaml_quote "${ARTIFACT_ROOT}")

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
      file_id_map: $(yaml_quote "${FIGURE_FILE_ID_MAP}")
      gene_map: $(yaml_quote "${FIGURE_GENE_MAP}")
      geneset_dir: $(yaml_quote "${FIGURE_GENESET_DIR}")
      gene_annotation: $(yaml_quote "${FIGURE_GENE_ANNOTATION}")
      posterior_dir: $(yaml_quote "${FIGURE_POSTERIOR_DIR}")
EOF
    if [[ -n "${FIGURE_POSTERIOR_NAME_MAP}" ]]; then
        cat >> "${config_path}" <<EOF
      posterior_name_map: $(yaml_quote "${FIGURE_POSTERIOR_NAME_MAP}")
EOF
    fi
}

# ============================================================
# 预计算基因集 YAML 列表（只算一次，避免每个 trait 重复 tr/printf 循环）
# ============================================================
GENESETS_YAML=""
IFS=',' read -r -a _gs_arr <<< "${FIGURE_GENESET_LIST}"
for _gs in "${_gs_arr[@]}"; do
    [[ -z "${_gs}" ]] && continue
    _gs="${_gs//\'/\'\'}"
    GENESETS_YAML+="          - '${_gs}'"$'\n'
done

# ============================================================
# 每种图类型的 YAML 配置片段生成
# ============================================================

append_cnmf_section() {
    local config_path="$1" lof_id="$2" trait_file="$3"
    cat >> "${config_path}" <<EOF
    cnmf:
      inputs:
        program_association_dir: $(yaml_quote "${FIGURE_PROGRAM_ASSOCIATION_DIR}")
        cnmf_regulation_dir: $(yaml_quote "${FIGURE_CNMF_REGULATION_DIR}")
      parameters:
        k: ${FIGURE_CNMF_K}
        trait_targets:
          - trait_file: $(yaml_quote "${trait_file}")
            trait_id: $(yaml_quote "${lof_id}")
EOF
}

append_burden_volcano_section() {
    local config_path="$1" lof_id="$2"
    cat >> "${config_path}" <<YAML_EOF
    burden_volcano:
      parameters:
        lof_ids: [$(yaml_quote "${lof_id}")]
        highlight_genesets:
${GENESETS_YAML}
YAML_EOF
}

append_gene_level_scatter_section() {
    local config_path="$1" lof_id="$2"
    cat >> "${config_path}" <<EOF
    gene_level_scatter:
      inputs:
        limma_path: $(yaml_quote "${FIGURE_LIMMA_PATH}")
        shet_path: $(yaml_quote "${FIGURE_SHET_PATH}")
      parameters:
        lof_ids: [$(yaml_quote "${lof_id}")]
        top_n_labels: ${FIGURE_GENE_LEVEL_TOP_N}
        y_limit: ${FIGURE_GENE_LEVEL_Y_LIMIT}
EOF
}

append_gene_level_qq_section() {
    local config_path="$1" lof_id="$2"
    cat >> "${config_path}" <<EOF
    gene_level_qq:
      inputs:
        limma_path: $(yaml_quote "${FIGURE_LIMMA_PATH}")
        shet_path: $(yaml_quote "${FIGURE_SHET_PATH}")
      parameters:
        lof_ids: [$(yaml_quote "${lof_id}")]
        y_limit: ${FIGURE_QQ_Y_LIMIT}
EOF
}

append_program_rankings_section() {
    local config_path="$1" lof_id="$2" trait_file="$3"
    cat >> "${config_path}" <<EOF
    program_rankings:
      inputs:
        program_association_dir: $(yaml_quote "${FIGURE_PROGRAM_ASSOCIATION_DIR}")
      parameters:
        k: ${FIGURE_CNMF_K}
        top_n: ${FIGURE_PROGRAM_TOP_N}
        trait_targets:
          - trait_file: $(yaml_quote "${trait_file}")
            trait_id: $(yaml_quote "${lof_id}")
EOF
}

append_gwas_manhattan_section() {
    local config_path="$1" gwas_id="$2"
    cat >> "${config_path}" <<YAML_EOF
    gwas_manhattan:
      parameters:
        gwas_ids: [$(yaml_quote "${gwas_id}")]
        flank_bp: ${FIGURE_GWAS_FLANK_BP}
        label_p_threshold: ${FIGURE_GWAS_LABEL_P_THRESHOLD}
        genomewide_threshold: ${FIGURE_GWAS_GENOMEWIDE_THRESHOLD}
        highlight_genesets:
${GENESETS_YAML}
YAML_EOF
}

# ============================================================
# Worker 函数：脚本自调用模式，避免 export -f 传递大量函数
# ============================================================
_fig_gen_one_lof() {
    local lof_id="$1"
    local trait_file lof_dir config_path script_path
    trait_file="$(resolve_trait_file "${lof_id}")"
    lof_dir="${LOF_ROOT}/${lof_id}"
    mkdir -p "${lof_dir}"
    config_path="${lof_dir}/${CONFIG_NAME}"

    write_config_header "${config_path}"
    if is_truthy "${FIGURE_ENABLE_CNMF}"; then
        append_cnmf_section "${config_path}" "${lof_id}" "${trait_file}"
    fi
    if is_truthy "${FIGURE_ENABLE_BURDEN_VOLCANO}"; then
        append_burden_volcano_section "${config_path}" "${lof_id}"
    fi
    if is_truthy "${FIGURE_ENABLE_GENE_LEVEL_SCATTER}"; then
        append_gene_level_scatter_section "${config_path}" "${lof_id}"
    fi
    if is_truthy "${FIGURE_ENABLE_PROGRAM_RANKINGS}"; then
        append_program_rankings_section "${config_path}" "${lof_id}" "${trait_file}"
    fi
    if is_truthy "${FIGURE_ENABLE_GENE_LEVEL_QQ}"; then
        append_gene_level_qq_section "${config_path}" "${lof_id}"
    fi

    if is_truthy "${FIGURE_ENABLE_CNMF}"; then
        script_path="${lof_dir}/figures_cnmf_${lof_id}.sh"
        write_slurm_script "${script_path}" "fig_cnmf_${lof_id}" \
            "${FIGURE_CNMF_MEM}" "${FIGURE_CNMF_CPUS}" \
            "${FIGURE_CNMF_PARTITION}" "${FIGURE_CNMF_TIME}" \
            "figures-cnmf" "${config_path}" "cnmf"
    fi
    if is_truthy "${FIGURE_ENABLE_BURDEN_VOLCANO}"; then
        script_path="${lof_dir}/figures_bvolcano_${lof_id}.sh"
        write_slurm_script "${script_path}" "fig_bvolcano_${lof_id}" \
            "${FIGURE_BURDEN_VOLCANO_MEM}" "${FIGURE_BURDEN_VOLCANO_CPUS}" \
            "${FIGURE_BURDEN_VOLCANO_PARTITION}" "${FIGURE_BURDEN_VOLCANO_TIME}" \
            "figures-burden-volcano" "${config_path}" "burden_volcano"
    fi
    if is_truthy "${FIGURE_ENABLE_GENE_LEVEL_SCATTER}"; then
        script_path="${lof_dir}/figures_gscatter_${lof_id}.sh"
        write_slurm_script "${script_path}" "fig_gscatter_${lof_id}" \
            "${FIGURE_GENE_LEVEL_SCATTER_MEM}" "${FIGURE_GENE_LEVEL_SCATTER_CPUS}" \
            "${FIGURE_GENE_LEVEL_SCATTER_PARTITION}" "${FIGURE_GENE_LEVEL_SCATTER_TIME}" \
            "figures-gene-level-scatter" "${config_path}" "gene_level_scatter"
    fi
    if is_truthy "${FIGURE_ENABLE_PROGRAM_RANKINGS}"; then
        script_path="${lof_dir}/figures_prank_${lof_id}.sh"
        write_slurm_script "${script_path}" "fig_prank_${lof_id}" \
            "${FIGURE_PROGRAM_RANKINGS_MEM}" "${FIGURE_PROGRAM_RANKINGS_CPUS}" \
            "${FIGURE_PROGRAM_RANKINGS_PARTITION}" "${FIGURE_PROGRAM_RANKINGS_TIME}" \
            "figures-program-rankings" "${config_path}" "program_rankings"
    fi
    if is_truthy "${FIGURE_ENABLE_GENE_LEVEL_QQ}"; then
        script_path="${lof_dir}/figures_glqq_${lof_id}.sh"
        write_slurm_script "${script_path}" "fig_glqq_${lof_id}" \
            "${FIGURE_GENE_LEVEL_QQ_MEM}" "${FIGURE_GENE_LEVEL_QQ_CPUS}" \
            "${FIGURE_GENE_LEVEL_QQ_PARTITION}" "${FIGURE_GENE_LEVEL_QQ_TIME}" \
            "figures-gene-level-qq" "${config_path}" "gene_level_qq"
    fi
}

_fig_gen_one_gwas() {
    local gwas_id="$1"
    local gwas_dir config_path script_path
    gwas_dir="${GWAS_ROOT}/${gwas_id}"
    mkdir -p "${gwas_dir}"
    config_path="${gwas_dir}/${CONFIG_NAME}"

    if is_truthy "${FIGURE_ENABLE_GWAS_MANHATTAN}"; then
        write_config_header "${config_path}"
        append_gwas_manhattan_section "${config_path}" "${gwas_id}"
        script_path="${gwas_dir}/figures_gman_${gwas_id}.sh"
        write_slurm_script "${script_path}" "fig_gman_${gwas_id}" \
            "${FIGURE_GWAS_MANHATTAN_MEM}" "${FIGURE_GWAS_MANHATTAN_CPUS}" \
            "${FIGURE_GWAS_MANHATTAN_PARTITION}" "${FIGURE_GWAS_MANHATTAN_TIME}" \
            "figures-gwas-manhattan" "${config_path}" "gwas_manhattan"
    fi
}

# Worker 入口（脚本自调用时通过这里进入）
case "${1:-}" in
    --worker-lof) _fig_gen_one_lof "$2"; exit 0 ;;
    --worker-gwas) _fig_gen_one_gwas "$2"; exit 0 ;;
esac

# ============================================================
# 清理并创建输出目录
# ============================================================
for _dir in "${LOF_ROOT}" "${GWAS_ROOT}" "${GROUP_ROOT}" "${PAIR_ROOT}" "${LOCUS_ROOT}"; do
    if [[ -z "${_dir}" || "${_dir}" == "/" ]]; then
        echo "FATAL: refusing to rm -rf empty or root path" >&2; exit 1
    fi
done
rm -rf "${LOF_ROOT}" "${GWAS_ROOT}" "${GROUP_ROOT}" "${PAIR_ROOT}" "${LOCUS_ROOT}"
mkdir -p "${LOF_ROOT}" "${GWAS_ROOT}" "${GROUP_ROOT}" "${PAIR_ROOT}" "${LOCUS_ROOT}"

# ============================================================
# 收集 ID 列表
# ============================================================
lof_ids=()
if [[ -n "${FIGURE_LOF_IDS}" ]]; then
    array_from_csv "${FIGURE_LOF_IDS}" lof_ids
else
    while IFS= read -r value; do
        lof_ids+=("${value}")
    done < <(collect_map_ids 2)
fi
if (( ${#lof_ids[@]} == 0 )); then
    echo "No LoF IDs were resolved from FIGURE_FILE_ID_MAP" >&2
    exit 1
fi

gwas_ids=()
if [[ -n "${FIGURE_GWAS_IDS}" ]]; then
    array_from_csv "${FIGURE_GWAS_IDS}" gwas_ids
else
    while IFS= read -r value; do
        gwas_ids+=("${value}")
    done < <(collect_map_ids 1)
fi

# ============================================================
# 提交队列
# ============================================================
declare -a submit_lines=()
declare -a lof_submit_lines=()
declare -a gwas_submit_lines=()
declare -a group_submit_lines=()
declare -a pair_submit_lines=()
declare -a locus_submit_lines=()

# ============================================================
# 1. 每个 LoF ID 的图（xargs -P 并行生成）
# ============================================================
echo "Generating LoF figure scripts (parallel ${FIGURE_PARALLEL})..."
printf '%s\n' "${lof_ids[@]}" | xargs -P "${FIGURE_PARALLEL}" -I {} \
    bash "$0" --worker-lof {}
echo "Done. Scanning LoF sbatch scripts..."
while IFS= read -r -d '' script_path; do
    submit_lines+=("sbatch \"${script_path}\"")
    lof_submit_lines+=("sbatch \"${script_path}\"")
done < <(find "${LOF_ROOT}" -name "figures_*.sh" -print0 | sort -z)

# ============================================================
# 2. 每个 GWAS ID 的图（xargs -P 并行生成）
# ============================================================
if (( ${#gwas_ids[@]} > 0 )); then
    echo "Generating GWAS figure scripts (parallel ${FIGURE_PARALLEL})..."
    printf '%s\n' "${gwas_ids[@]}" | xargs -P "${FIGURE_PARALLEL}" -I {} \
        bash "$0" --worker-gwas {}
    echo "Done. Scanning GWAS sbatch scripts..."
    while IFS= read -r -d '' script_path; do
        submit_lines+=("sbatch \"${script_path}\"")
        gwas_submit_lines+=("sbatch \"${script_path}\"")
    done < <(find "${GWAS_ROOT}" -name "figures_*.sh" -print0 | sort -z)
fi

# ============================================================
# 3. 汇总图（跨 trait）
# ============================================================
overview_ids=()
if [[ -n "${FIGURE_HEATMAP_LOF_IDS}" ]]; then
    array_from_csv "${FIGURE_HEATMAP_LOF_IDS}" overview_ids
else
    overview_ids=("${lof_ids[@]}")
fi

# --- Program Heatmap ---
if is_truthy "${FIGURE_ENABLE_PROGRAM_HEATMAP}" && (( ${#overview_ids[@]} > 0 )); then
    heatmap_dir="${GROUP_ROOT}/program_heatmap"
    mkdir -p "${heatmap_dir}"
    config_path="${heatmap_dir}/${CONFIG_NAME}"
    write_config_header "${config_path}"
    heatmap_id="$(batch_id_from_ids "program_heatmap" "${overview_ids[@]}")"
    cat >> "${config_path}" <<EOF
    program_heatmap:
      inputs:
        program_association_dir: $(yaml_quote "${FIGURE_PROGRAM_ASSOCIATION_DIR}")
      parameters:
        k: ${FIGURE_CNMF_K}
        output_id: $(yaml_quote "${heatmap_id}")
        trait_targets:
EOF
    for lof_id in "${overview_ids[@]}"; do
        trait_file="$(resolve_trait_file "${lof_id}")"
        cat >> "${config_path}" <<EOF
          - trait_file: $(yaml_quote "${trait_file}")
            trait_id: $(yaml_quote "${lof_id}")
EOF
    done
    script_path="${heatmap_dir}/figures_program_heatmap.sh"
    write_slurm_script "${script_path}" "fig_heatmap" \
        "${FIGURE_PROGRAM_HEATMAP_MEM}" "${FIGURE_PROGRAM_HEATMAP_CPUS}" \
        "${FIGURE_PROGRAM_HEATMAP_PARTITION}" "${FIGURE_PROGRAM_HEATMAP_TIME}" \
        "figures" "${config_path}" "program_heatmap"
    submit_lines+=("sbatch \"${script_path}\"")
    group_submit_lines+=("sbatch \"${script_path}\"")
fi


# --- Cross-Trait Heatmap ---
if is_truthy "${FIGURE_ENABLE_CROSS_TRAIT_HEATMAP}" && (( ${#overview_ids[@]} > 1 )); then
    cth_dir="${GROUP_ROOT}/cross_trait_heatmap"
    mkdir -p "${cth_dir}"
    config_path="${cth_dir}/${CONFIG_NAME}"
    write_config_header "${config_path}"
    cross_heatmap_id="$(batch_id_from_ids "cross_trait_heatmap" "${overview_ids[@]}")"
    cat >> "${config_path}" <<EOF
    cross_trait_heatmap:
      parameters:
        output_id: $(yaml_quote "${cross_heatmap_id}")
        method: $(yaml_quote "${FIGURE_CROSS_TRAIT_METHOD}")
        lof_ids:
EOF
    for lof_id in "${overview_ids[@]}"; do
        cat >> "${config_path}" <<EOF
          - $(yaml_quote "${lof_id}")
EOF
    done
    script_path="${cth_dir}/figures_cross_trait_heatmap.sh"
    write_slurm_script "${script_path}" "fig_cth" \
        "${FIGURE_CROSS_TRAIT_HEATMAP_MEM}" "${FIGURE_CROSS_TRAIT_HEATMAP_CPUS}" \
        "${FIGURE_CROSS_TRAIT_HEATMAP_PARTITION}" "${FIGURE_CROSS_TRAIT_HEATMAP_TIME}" \
        "figures" "${config_path}" "cross_trait_heatmap"
    submit_lines+=("sbatch \"${script_path}\"")
    group_submit_lines+=("sbatch \"${script_path}\"")
fi

# --- cNMF Program Top Genes ---
if is_truthy "${FIGURE_ENABLE_CNMF_PROGRAM_TOP_GENES}"; then
    topg_dir="${GROUP_ROOT}/cnmf_program_top_genes"
    mkdir -p "${topg_dir}"
    config_path="${topg_dir}/${CONFIG_NAME}"
    write_config_header "${config_path}"
    cat >> "${config_path}" <<EOF
    cnmf_program_top_genes:
      inputs:
        spectra_path: $(yaml_quote "${FIGURE_SPECTRA_PATH}")
      parameters:
        k: ${FIGURE_CNMF_K}
EOF
    script_path="${topg_dir}/figures_cnmf_top_genes.sh"
    write_slurm_script "${script_path}" "fig_topg" \
        "${FIGURE_CNMF_PROGRAM_TOP_GENES_MEM}" "${FIGURE_CNMF_PROGRAM_TOP_GENES_CPUS}" \
        "${FIGURE_CNMF_PROGRAM_TOP_GENES_PARTITION}" "${FIGURE_CNMF_PROGRAM_TOP_GENES_TIME}" \
        "figures" "${config_path}" "cnmf_program_top_genes"
    submit_lines+=("sbatch \"${script_path}\"")
    group_submit_lines+=("sbatch \"${script_path}\"")
fi

# --- cNMF Program Enrichment (气泡图, 54 基因集) ---
if is_truthy "${FIGURE_ENABLE_CNMF_PROGRAM_ENRICHMENT}"; then
    enrich_dir="${GROUP_ROOT}/cnmf_program_enrichment"
    mkdir -p "${enrich_dir}"
    config_path="${enrich_dir}/${CONFIG_NAME}"
    write_config_header "${config_path}"
    cat >> "${config_path}" <<YAML_EOF
    cnmf_program_enrichment:
      inputs:
        spectra_path: $(yaml_quote "${FIGURE_SPECTRA_PATH}")
      parameters:
        k: ${FIGURE_CNMF_K}
        genesets:
${GENESETS_YAML}
YAML_EOF
    script_path="${enrich_dir}/figures_cnmf_enrichment.sh"
    write_slurm_script "${script_path}" "fig_enrich" \
        "${FIGURE_CNMF_PROGRAM_ENRICHMENT_MEM}" "${FIGURE_CNMF_PROGRAM_ENRICHMENT_CPUS}" \
        "${FIGURE_CNMF_PROGRAM_ENRICHMENT_PARTITION}" "${FIGURE_CNMF_PROGRAM_ENRICHMENT_TIME}" \
        "figures" "${config_path}" "cnmf_program_enrichment"
    submit_lines+=("sbatch \"${script_path}\"")
    group_submit_lines+=("sbatch \"${script_path}\"")
fi

# ============================================================
# 4. Cross-Trait 成对图
# ============================================================
if is_truthy "${FIGURE_ENABLE_CROSS_TRAIT}" && [[ -n "${FIGURE_CROSS_TRAIT_PAIRS}" ]]; then
    IFS=',' read -r -a cross_pairs <<< "${FIGURE_CROSS_TRAIT_PAIRS}"
    for raw_pair in "${cross_pairs[@]}"; do
        raw_pair="$(trim "${raw_pair}")"
        [[ -z "${raw_pair}" ]] && continue
        IFS=':' read -r id_x id_y <<< "${raw_pair}"
        id_x="$(trim "${id_x:-}")"
        id_y="$(trim "${id_y:-}")"
        [[ -n "${id_x}" && -n "${id_y}" ]] || continue
        pair_id="$(safe_name "${id_x}__${id_y}")"
        pair_dir="${PAIR_ROOT}/${pair_id}"
        mkdir -p "${pair_dir}"
        config_path="${pair_dir}/${CONFIG_NAME}"
        write_config_header "${config_path}"
        cat >> "${config_path}" <<EOF
    cross_trait:
      parameters:
        top_n_labels: ${FIGURE_CROSS_TRAIT_TOP_N}
        lof_pairs:
          - id_x: $(yaml_quote "${id_x}")
            id_y: $(yaml_quote "${id_y}")
            output_id: $(yaml_quote "${pair_id}")
EOF
        script_path="${pair_dir}/figures_cross_trait_${pair_id}.sh"
        write_slurm_script "${script_path}" "figx_${pair_id}" \
            "${FIGURE_CROSS_TRAIT_MEM}" "${FIGURE_CROSS_TRAIT_CPUS}" \
            "${FIGURE_CROSS_TRAIT_PARTITION}" "${FIGURE_CROSS_TRAIT_TIME}" \
            "figures" "${config_path}" "cross_trait"
        submit_lines+=("sbatch \"${script_path}\"")
        pair_submit_lines+=("sbatch \"${script_path}\"")
    done
fi

# ============================================================
# 5. GWAS Locus Zoom
# ============================================================
run_locus_zoom=0
case "${FIGURE_ENABLE_GWAS_LOCUS_ZOOM,,}" in
    auto|'')
        [[ -n "${FIGURE_LOCI_TSV}" ]] && run_locus_zoom=1
        ;;
    1|true|yes|on)
        run_locus_zoom=1
        ;;
    0|false|no|off)
        run_locus_zoom=0
        ;;
    *)
        echo "FIGURE_ENABLE_GWAS_LOCUS_ZOOM must be auto or a boolean value. Got: ${FIGURE_ENABLE_GWAS_LOCUS_ZOOM}" >&2
        exit 1
        ;;
esac

if (( run_locus_zoom == 1 )); then
    [[ -n "${FIGURE_LOCI_TSV}" ]] || { echo "FIGURE_LOCI_TSV is required when FIGURE_ENABLE_GWAS_LOCUS_ZOOM is enabled" >&2; exit 1; }
    FIGURE_LOCI_TSV="$(resolve_path "${FIGURE_LOCI_TSV}")"
    [[ -f "${FIGURE_LOCI_TSV}" ]] || { echo "FIGURE_LOCI_TSV does not exist: ${FIGURE_LOCI_TSV}" >&2; exit 1; }

    while IFS=$'\x1f' read -r gwas_id gwas_file chrom start end locus_label flank_bp output_id; do
        gwas_id="$(trim "${gwas_id:-}")"
        gwas_file="$(trim "${gwas_file:-}")"
        [[ -n "${gwas_id}" || -n "${gwas_file}" ]] || continue
        chrom="$(trim "${chrom:-}")"
        start="$(trim "${start:-}")"
        end="$(trim "${end:-}")"
        [[ -n "${chrom}" && -n "${start}" && -n "${end}" ]] || continue
        locus_label="${locus_label:-chr${chrom}:${start}-${end}}"
        flank_bp="${flank_bp:-${FIGURE_LOCUS_FLANK_BP}}"
        if [[ -z "${output_id:-}" ]]; then
            base_id="${gwas_id:-$(safe_name "$(basename "${gwas_file}")")}"
            output_id="${base_id}__chr${chrom}_${start}_${end}"
        fi
        output_id="$(safe_name "${output_id}")"
        locus_dir="${LOCUS_ROOT}/${output_id}"
        mkdir -p "${locus_dir}"
        config_path="${locus_dir}/${CONFIG_NAME}"
        write_config_header "${config_path}"
        if [[ -n "${gwas_id}" ]]; then
            cat >> "${config_path}" <<EOF
    gwas_locus_zoom:
      parameters:
        genomewide_threshold: ${FIGURE_GWAS_GENOMEWIDE_THRESHOLD}
        label_top_n: ${FIGURE_LOCUS_LABEL_TOP_N}
        locus_targets:
          - gwas_id: $(yaml_quote "${gwas_id}")
            chrom: $(yaml_quote "${chrom}")
            start: ${start}
            end: ${end}
            flank_bp: ${flank_bp}
            locus_label: $(yaml_quote "${locus_label}")
            output_id: $(yaml_quote "${output_id}")
EOF
        elif [[ -n "${gwas_file}" ]]; then
            cat >> "${config_path}" <<EOF
    gwas_locus_zoom:
      parameters:
        genomewide_threshold: ${FIGURE_GWAS_GENOMEWIDE_THRESHOLD}
        label_top_n: ${FIGURE_LOCUS_LABEL_TOP_N}
        locus_targets:
          - gwas_file: $(yaml_quote "${gwas_file}")
            chrom: $(yaml_quote "${chrom}")
            start: ${start}
            end: ${end}
            flank_bp: ${flank_bp}
            locus_label: $(yaml_quote "${locus_label}")
            output_id: $(yaml_quote "${output_id}")
EOF
        else
            echo "Each FIGURE_LOCI_TSV row must provide gwas_id or gwas_file" >&2
            exit 1
        fi
        script_path="${locus_dir}/figures_gwas_locus_${output_id}.sh"
        write_slurm_script "${script_path}" "figl_${output_id}" \
            "${FIGURE_GWAS_LOCUS_ZOOM_MEM}" "${FIGURE_GWAS_LOCUS_ZOOM_CPUS}" \
            "${FIGURE_GWAS_LOCUS_ZOOM_PARTITION}" "${FIGURE_GWAS_LOCUS_ZOOM_TIME}" \
            "figures" "${config_path}" "gwas_locus_zoom"
        submit_lines+=("sbatch \"${script_path}\"")
        locus_submit_lines+=("sbatch \"${script_path}\"")
    done < <(
        awk -F'\t' 'NR > 1 {
            printf "%s\037%s\037%s\037%s\037%s\037%s\037%s\037%s\n", $1, $2, $3, $4, $5, $6, $7, $8
        }' "${FIGURE_LOCI_TSV}"
    )
fi

# ============================================================
# 生成提交脚本
# ============================================================
submit_script="${BATCH_ROOT}/submit_figures.sh"
submit_lof_script="${BATCH_ROOT}/submit_figures_lof.sh"
submit_gwas_script="${BATCH_ROOT}/submit_figures_gwas.sh"
submit_groups_script="${BATCH_ROOT}/submit_figures_groups.sh"
submit_pairs_script="${BATCH_ROOT}/submit_figures_pairs.sh"
submit_loci_script="${BATCH_ROOT}/submit_figures_loci.sh"
submit_aggregate_script="${BATCH_ROOT}/submit_figures_aggregate.sh"

write_submit_helper_from_arrays "${submit_script}" submit_lines
write_submit_helper_from_arrays "${submit_lof_script}" lof_submit_lines
write_submit_helper_from_arrays "${submit_gwas_script}" gwas_submit_lines
write_submit_helper_from_arrays "${submit_groups_script}" group_submit_lines
write_submit_helper_from_arrays "${submit_pairs_script}" pair_submit_lines
write_submit_helper_from_arrays "${submit_loci_script}" locus_submit_lines
write_submit_helper_from_arrays "${submit_aggregate_script}" \
    group_submit_lines \
    pair_submit_lines \
    locus_submit_lines

echo "Generated figure sbatch scripts under: ${BATCH_ROOT}"
echo "Submit helper: ${submit_script}"
