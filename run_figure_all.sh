#!/usr/bin/env bash
#
# figure_all 专用启动脚本
# 用法: bash run_figure_all.sh
# 放在 /gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/ 下运行
#
set -euo pipefail

# ============================================================
# 数据路径（三个来源）
# ============================================================
export PROJECT_ROOT=/gpfs/chencao/qinminzhang/Nature/mine/code
export ARTIFACT_ROOT=/gpfs/chencao/qinminzhang/Nature/mine/code/outputs

# file_id_map（4 列 TSV：id1, id2, path1, path2）
export FIGURE_FILE_ID_MAP=${PROJECT_ROOT}/configs/path.file_id_map.tsv

# posterior 目录（分 trait 跑的）
export FIGURE_POSTERIOR_DIR=/gpfs/chencao/qinminzhang/workflow/catalog_lof/run_all/outputs/genebayes/posterior

# cNMF 数据（跑一次的）
# cNMF regulation（跑一次）和 association（分 trait）在不同位置
export FIGURE_CNMF_REGULATION_DIR=/gpfs/chencao/qinminzhang/Nature/mine/code/outputs/perturbseq/cnmf_genomewide/cNMF_regulation/K562GW
export FIGURE_PROGRAM_ASSOCIATION_DIR=/gpfs/chencao/qinminzhang/workflow/catalog_lof/run_all/outputs/perturbseq/cnmf_genomewide/trait_association/K562GW/ProgramLevel

# limma 矩阵 + shet（gene_level_scatter 自己算 correlation，不依赖 plot）
export FIGURE_LIMMA_PATH=/gpfs/chencao/qinminzhang/Nature/mine/code/outputs/perturbseq/gene_level/K562GW/limma_logFC_sum.txt
export FIGURE_SHET_PATH=/gpfs/chencao/qinminzhang/Nature/mine/code/data/shet_10bins.txt

# 输出目录
export BATCH_ROOT=$(pwd)/scripts
export LOGS_ROOT=$(pwd)/logs/figures
export ARTIFACT_ROOT=$(pwd)/outputs

# conda 环境
export CONTROL_ENV=paper-pipeline-control
export FIGURES_R_ENV=paper-pipeline-plot
export CONDA_SH=${HOME}/miniconda3/etc/profile.d/conda.sh

# ============================================================
# 开关：全部启用
# ============================================================
export FIGURE_ENABLE_CNMF=1
export FIGURE_ENABLE_PROGRAM_RANKINGS=1
export FIGURE_ENABLE_PROGRAM_HEATMAP=1
export FIGURE_ENABLE_BURDEN_VOLCANO=1
export FIGURE_ENABLE_GWAS_MANHATTAN=1
export FIGURE_ENABLE_GENE_LEVEL_SCATTER=1
export FIGURE_ENABLE_GENE_LEVEL_QQ=1
export FIGURE_ENABLE_CROSS_TRAIT=0                         # pairs 需要手动配
export FIGURE_ENABLE_CROSS_TRAIT_HEATMAP=1
export FIGURE_ENABLE_GWAS_LOCUS_ZOOM=0                     # 需要 loci TSV
export FIGURE_ENABLE_CNMF_PROGRAM_TOP_GENES=1
export FIGURE_ENABLE_CNMF_PROGRAM_ENRICHMENT=1

# ============================================================
#  资源：每种图独立配置
#
#  可用分区:
#    cu,privority,batch01 : 56 cores, 256 GB RAM
#    fat                  : 56 cores,   1 TB RAM
#
#  保守策略: 所有图给到 cu,privority,batch01 (256G 足够),
#           只有 GWAS Manhattan 上 fat (1T) 防止大文件 OOM
# ============================================================
export FIGURE_DEFAULT_PARTITION=cu,privority,batch01

# --- cNMF 散点图 (60 programs, ~2K perturb genes) ---
export FIGURE_CNMF_MEM=8G
export FIGURE_CNMF_CPUS=2
export FIGURE_CNMF_PARTITION=cu,privority,batch01
export FIGURE_CNMF_TIME=12:00:00

# --- Burden 火山图 (~20K genes × 54 genesets) ---
export FIGURE_BURDEN_VOLCANO_MEM=8G
export FIGURE_BURDEN_VOLCANO_CPUS=2
export FIGURE_BURDEN_VOLCANO_PARTITION=cu,privority,batch01
export FIGURE_BURDEN_VOLCANO_TIME=12:00:00

# --- Gene-level 散点图 (posterior + correlation, ~20K genes) ---
export FIGURE_GENE_LEVEL_SCATTER_MEM=8G
export FIGURE_GENE_LEVEL_SCATTER_CPUS=2
export FIGURE_GENE_LEVEL_SCATTER_PARTITION=cu,privority,batch01
export FIGURE_GENE_LEVEL_SCATTER_TIME=12:00:00

# --- Program 排名条形图 (60 programs 排序, 最轻) ---
export FIGURE_PROGRAM_RANKINGS_MEM=4G
export FIGURE_PROGRAM_RANKINGS_CPUS=1
export FIGURE_PROGRAM_RANKINGS_PARTITION=cu,privority,batch01
export FIGURE_PROGRAM_RANKINGS_TIME=12:00:00

# --- GWAS Manhattan（★ 唯一上 fat 的, 受百万级变异位点 + 54 基因集注释 + 光栅化）---
export FIGURE_GWAS_MANHATTAN_MEM=128G
export FIGURE_GWAS_MANHATTAN_CPUS=8
export FIGURE_GWAS_MANHATTAN_PARTITION=fat
export FIGURE_GWAS_MANHATTAN_TIME=48:00:00

# --- Program × Trait 热图 (N trait × 60 programs + hclust) ---
export FIGURE_PROGRAM_HEATMAP_MEM=32G
export FIGURE_PROGRAM_HEATMAP_CPUS=4
export FIGURE_PROGRAM_HEATMAP_PARTITION=cu,privority,batch01
export FIGURE_PROGRAM_HEATMAP_TIME=24:00:00

# --- Gene-level QQ 图 (N trait × ~20K genes correlation) ---
export FIGURE_GENE_LEVEL_QQ_MEM=32G
export FIGURE_GENE_LEVEL_QQ_CPUS=4
export FIGURE_GENE_LEVEL_QQ_PARTITION=cu,privority,batch01
export FIGURE_GENE_LEVEL_QQ_TIME=24:00:00

# --- Cross-Trait 相关性热图 (N trait posterior merge + N×N cor) ---
export FIGURE_CROSS_TRAIT_HEATMAP_MEM=32G
export FIGURE_CROSS_TRAIT_HEATMAP_CPUS=4
export FIGURE_CROSS_TRAIT_HEATMAP_PARTITION=cu,privority,batch01
export FIGURE_CROSS_TRAIT_HEATMAP_TIME=24:00:00

# --- cNMF program top 基因条形图 ---
export FIGURE_CNMF_PROGRAM_TOP_GENES_MEM=8G
export FIGURE_CNMF_PROGRAM_TOP_GENES_CPUS=2
export FIGURE_CNMF_PROGRAM_TOP_GENES_PARTITION=cu,privority,batch01
export FIGURE_CNMF_PROGRAM_TOP_GENES_TIME=12:00:00

# --- cNMF 基因集富集气泡图 (60 programs × 54 genesets hypergeometric) ---
export FIGURE_CNMF_PROGRAM_ENRICHMENT_MEM=8G
export FIGURE_CNMF_PROGRAM_ENRICHMENT_CPUS=2
export FIGURE_CNMF_PROGRAM_ENRICHMENT_PARTITION=cu,privority,batch01
export FIGURE_CNMF_PROGRAM_ENRICHMENT_TIME=12:00:00

# --- Cross-Trait 散点图 (2 trait posterior, 默认关闭) ---
export FIGURE_CROSS_TRAIT_MEM=8G
export FIGURE_CROSS_TRAIT_CPUS=2
export FIGURE_CROSS_TRAIT_PARTITION=cu,privority,batch01
export FIGURE_CROSS_TRAIT_TIME=12:00:00

# --- GWAS Locus Zoom (单染色体区域, 默认关闭) ---
export FIGURE_GWAS_LOCUS_ZOOM_MEM=8G
export FIGURE_GWAS_LOCUS_ZOOM_CPUS=2
export FIGURE_GWAS_LOCUS_ZOOM_PARTITION=cu,privority,batch01
export FIGURE_GWAS_LOCUS_ZOOM_TIME=12:00:00

# ============================================================
# 执行生成器
# ============================================================
GENERATOR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/generate_figure_sbatch.sh

if [[ ! -f "${GENERATOR}" ]]; then
    echo "找不到 generate_figure_sbatch.sh" >&2
    echo "路径: ${GENERATOR}" >&2
    echo "请确认 generate_figure_sbatch.sh 已放到该路径" >&2
    exit 1
fi

echo "=== Figure 数据路径 ==="
echo "PROJECT_ROOT                = ${PROJECT_ROOT}"
echo "ARTIFACT_ROOT               = ${ARTIFACT_ROOT}"
echo "FIGURE_FILE_ID_MAP          = ${FIGURE_FILE_ID_MAP}"
echo "FIGURE_POSTERIOR_DIR        = ${FIGURE_POSTERIOR_DIR}"
echo "FIGURE_PROGRAM_ASSOC_DIR    = ${FIGURE_PROGRAM_ASSOCIATION_DIR}"
echo "FIGURE_CNMF_REGULATION_DIR  = ${FIGURE_CNMF_REGULATION_DIR}"
echo "FIGURE_LIMMA_PATH          = ${FIGURE_LIMMA_PATH}"
echo "FIGURE_SHET_PATH           = ${FIGURE_SHET_PATH}"
echo "BATCH_ROOT                  = ${BATCH_ROOT}"
echo "LOGS_ROOT                   = ${LOGS_ROOT}"
echo ""

bash "${GENERATOR}"
