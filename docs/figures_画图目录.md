# Figures 画图目录

本文档列出当前所有可画的图像，每张图包含：图像描述、输入数据来源、核心画图代码。代码来自 `src/paper_pipeline/figures/scripts/` 下的 R 脚本。

画图风格统一为 `theme_classic(base_family = "Helvetica")`，颜色调色板为 `#B40426`(红) / `#3B4CC0`(蓝) / `#34A853`(绿) / `#FEA601`(橙)。

---

## 一、单 Trait 图（一个 trait 一张图）

### 1. Program 关联散点图 (cnmf_program_regulator_scatter)

**数据来源：** cNMF association 目录下的两个文件

| 输入文件 | 必需字段 |
|----------|----------|
| `programs_enrichment_K{k}_{trait_file}` | `Program`, `MEANgamma_top100`, `shet_adjusted_random_mean`, `MEANgamma_top100_shet_adjusted_P` |
| `regulators_enrichment_K{k}_{trait_file}` | `Program`, `beta_withShet`, `P_withShet` |

**数据处理（helpers.R `read_program_regulator_summary`）：**
```r
# Merge 两个 enrichment 表，计算指标
bonferroni_cutoff <- 0.05 / nrow(df)
df$program_score <- sign(MEANgamma_top100 - shet_adjusted_random_mean) * (-log10(MEANgamma_top100_shet_adjusted_P))
df$regulator_score <- sign(beta_withShet) * (-log10(P_withShet))
# 着色：both_enriched / program_enriched / regulator_enriched / other
```

**核心画图代码：**
```r
library(ggplot2); library(ggrepel)

g <- ggplot(df, aes(x = program_score, y = regulator_score, label = label, color = color))
g <- g + theme_classic(base_size = 24, base_family = "Helvetica")
g <- g + geom_point(size = 4)
g <- g + geom_vline(xintercept = 0, linetype = "dashed")
g <- g + geom_hline(yintercept = 0, linetype = "dashed")
g <- g + geom_text_repel(size = 7, max.overlaps = 100)
g <- g + scale_color_manual(values = c(
  "other" = "grey60",
  "program_enriched" = "#FEA601",
  "regulator_enriched" = "#4783B5",
  "both_enriched" = "#34A853"
))
g <- g + xlab("Program burden effect, signed -log10(P)")
g <- g + ylab("Regulator-burden correlation, signed -log10(P)")
```

---

### 2. LoF Burden 火山图 (burden_volcano)

**数据来源：**

| 输入文件 | 必需字段 |
|----------|----------|
| LoF burden CSV (如 `*_M1_001.summary_statistics.csv`) | `beta`, `standard_error`, `ensg` |
| `gencode_v41_gname_gid_ALL_sorted_onlyID` | 2列无表头：`ENSG` `gene_symbol` |
| `geneset_dir/{geneset}.txt` | 每行一个 gene symbol |

**数据处理：**
```r
lof <- data.table::fread(burden_path, data.table = FALSE)
lof$P <- 2 * pnorm(-abs(lof$beta / lof$standard_error))
df <- lof[, c("beta", "P", "ensg")]
df$FDR <- p.adjust(df$P, method = "BH")
df$LABEL <- label_from_ensg(df$ensg, gene_map_path)
df$neg_log10_p <- -log10(pmax(df$P, .Machine$double.xmin))
# FDR 线：满足 FDR <= line_fdr_threshold 的最大 P 值
```

**核心画图代码：**
```r
library(ggplot2); library(ggrepel)

g <- ggplot(df, aes(x = beta, y = neg_log10_p, color = geneset))
g <- g + theme_classic(base_size = 20, base_family = "Helvetica")
g <- g + geom_point(alpha = 0.45, size = 2.5)
g <- g + geom_point(data = highlight_df, alpha = 0.9, size = 3)     # 基因集内点加亮
g <- g + geom_text_repel(data = label_df, aes(label = label),        # 标注 FDR 显著的基因集基因
    size = 5, max.overlaps = 100, show.legend = FALSE)
g <- g + geom_hline(yintercept = threshold_y, linetype = "dashed")   # FDR 线
g <- g + scale_color_manual(values = palette_values)                 # other=grey70, 基因集=蓝红绿橙
g <- g + xlab("Effect size")
g <- g + ylab("-log10(P)")
```

---

### 3. GWAS Manhattan 图 (gwas_manhattan)

**数据来源：**

| 输入文件 | 必需字段 |
|----------|----------|
| GWAS summary (如 `{trait}_irnt.tsv.gz`) | `variant` (CHROM:POS 格式), `pval`, 可选 `low_confidence_variant` |
| `genes.protein_coding.v39.gtf` | 标准 GTF 或 4列表：`gene`, `chr`, `start`, `end` |
| `geneset_dir/{geneset}.txt` | 每行一个 gene symbol |

**数据处理（helpers.R `parse_variant_table` + `add_cumulative_coordinates`）：**
```r
# 解析 variant = "1:55505221:A:G" → CHROM="1", POS=55505221
# 只保留 1-22 号染色体，过滤 low_confidence_variant==TRUE
# 计算累加坐标：COORD = POS + 前染色体最大坐标偏移
coord_data <- add_cumulative_coordinates(gwas_df)
# 生成染色体轴刻度：axis_df = data.frame(CHROM=c("1","2"...), center=c(中位坐标))
```

**核心画图代码：**
```r
library(ggplot2); library(ggrepel); library(ggrastr)

g <- ggplot(gwas_df, aes(x = COORD, y = neg_log10_p))
g <- g + theme_classic(base_size = 20, base_family = "Helvetica")
g <- g + ggrastr::geom_point_rast(color = "grey70", alpha = 0.5, size = 1.2)  # 光栅化背景
g <- g + geom_hline(yintercept = -log10(genomewide_threshold), linetype = "dashed", color = "red")
g <- g + geom_point(data = highlight_hits, aes(color = geneset), size = 2.5)  # 基因集周边高亮
g <- g + geom_text_repel(data = labeled_hits, aes(label = LABEL, color = geneset),
    size = 5, max.overlaps = 100)
g <- g + scale_x_continuous(labels = axis_df$CHROM, breaks = axis_df$center)  # 染色体刻度
g <- g + scale_color_manual(values = palette_values)                          # other=grey70
g <- g + xlab("Chromosome")
g <- g + ylab("-log10(P)")
```

---

### 4. Gene-Level 散点图 (gene_level_scatter)

**数据来源：**

| 输入文件 | 必需字段 |
|----------|----------|
| `{trait_stem}.per_gene_estimates.tsv` | `ensg`, `post_mean` |
| `{trait_stem}..._geneRegulation_correlation.txt` | `ensg`, `P_withShet`, `beta_withShet`，可选 `gene` |
| `gencode_v41_gname_gid_ALL_sorted_onlyID` | 2列无表头 |

**数据处理：**
```r
# merge 两张表 by ensg
df$signed_log10_p <- sign(beta_withShet) * (-log10(pmax(P_withShet, .Machine$double.xmin)))
df$fdr <- p.adjust(df$P_withShet, method = "BH")
# label：highlight_genes 中指定的 + top_n 按 |post_mean * signed_log10_p| 最大的基因
df$label_score <- abs(post_mean) * abs(signed_log10_p)
# FDR 0.05 线
```

**核心画图代码：**
```r
library(ggplot2); library(ggrepel)

g <- ggplot(df, aes(x = post_mean, y = signed_log10_p))
g <- g + theme_classic(base_size = 18, base_family = "Helvetica")
g <- g + geom_point(alpha = 0.6, size = 2.5, color = "#4C78A8")
g <- g + geom_hline(yintercept = 0, linetype = "dashed", color = "grey60")
g <- g + geom_vline(xintercept = 0, linetype = "dashed", color = "grey60")
g <- g + geom_hline(yintercept = c(threshold_y, -threshold_y), linetype = "dotted", color = "#B40426")
g <- g + geom_text_repel(data = label_df, aes(label = label),
    size = 4.5, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.2, min.segment.length = 0)
g <- g + xlab(sprintf("%s posterior effect", trait_label))
g <- g + ylab("Gene-level signed -log10(P)")
g <- g + coord_cartesian(ylim = c(-y_limit, y_limit))
```

---

### 5. Program 排名条形图 (program_rankings)

**数据来源：** 与 #1 相同（programs_enrichment + regulators_enrichment）

**核心画图代码：**
```r
library(ggplot2)

# 对 program_score 和 regulator_score 各画一张水平条形图
plot_top_metric <- function(metric_name, rank_name, title_text) {
  tmp <- head(df[order(df[[rank_name]]), ], top_n)
  tmp$Program <- factor(tmp$Program, levels = rev(tmp$Program))  # 反转使 top1 在顶部

  g <- ggplot(tmp, aes(x = Program, y = value, fill = value > 0))
  g <- g + theme_classic(base_size = 18, base_family = "Helvetica")
  g <- g + geom_col(width = 0.75)
  g <- g + coord_flip()                                           # 水平条形
  g <- g + geom_hline(yintercept = 0, linetype = "dashed")
  g <- g + scale_fill_manual(values = c("TRUE" = "#B40426", "FALSE" = "#3B4CC0"), guide = "none")
  g <- g + xlab("Program")
  g <- g + ylab(title_text)
}
# 调用两次：program_score 和 regulator_score
```

---

### 6. GWAS Locus Zoom (gwas_locus_zoom)

**数据来源：** 与 #3 相同（GWAS + gene annotation），额外指定 `chrom`, `start`, `end`, `flank_bp`

**数据处理：**
```r
# 按染色体和坐标区间截取 GWAS 数据
region_df <- gwas_df[gwas_df$CHROM == chrom & POS >= start - flank_bp & POS <= end + flank_bp, ]
# 从 gene annotation 截取区间内基因
gene_df <- annotation[annotation$chr == chrom & end >= region_start & start <= region_end, ]
# 基因轨道分层渲染，避免重叠：track_index = (row_index - 1) %% 3 → 3层
```

**核心画图代码：**
```r
library(ggplot2); library(ggrepel)

g <- ggplot(region_df, aes(x = POS, y = neg_log10_p))
g <- g + theme_classic(base_size = 18, base_family = "Helvetica")
g <- g + geom_point(aes(color = in_core_region), alpha = 0.7, size = 2.2)
g <- g + scale_color_manual(values = c("TRUE" = "#B40426", "FALSE" = "grey65"), guide = "none")
g <- g + geom_hline(yintercept = -log10(genomewide_threshold), linetype = "dashed", color = "#B40426")
g <- g + geom_vline(xintercept = c(start, end), linetype = "dotted", color = "grey50")
# 基因轨道
g <- g + geom_segment(data = gene_track_df, aes(x = start, xend = end, y = track_y, yend = track_y),
    linewidth = 1.1, color = "#2F4F4F")
g <- g + geom_text(data = gene_track_df, aes(x = (start+end)/2, y = track_y, label = gene),
    vjust = 1.6, size = 3.6, color = "#2F4F4F")
# 标注 top hits
g <- g + geom_text_repel(data = top_hits, aes(label = label), size = 4, max.overlaps = Inf)
g <- g + xlab(sprintf("chr%s position", chrom))
g <- g + ylab("-log10(P)")
```

---

### 7. trans-eQTL vs LoF 散点图 (legacy Fig4J)

**数据来源：**

| 输入文件 | 必需字段 |
|----------|----------|
| `regulators_enrichment_K{k}_{trait_file}` | `Program`, `beta_withShet`, `P_withShet` |
| `program_trans/result/P{i}_transeQTL_*_directionalTest.txt` | `mean_Z_GWAS`, `mean_Z_ctrl`, `Ttest_P` |

**核心画图代码：**
```r
# 每个 program：LOF score = signed -log10(P_withShet)，GWAS score = signed -log10(Ttest_P)
# 散点图 + ggrepel 标注
g <- ggplot(df, aes(x = LOFscore, y = GWASscore, label = LABEL))
g <- g + geom_point(size = 3)
g <- g + geom_vline(xintercept = 0, linetype = "dashed")
g <- g + geom_hline(yintercept = 0, linetype = "dashed")
g <- g + geom_text_repel(size = 5, max.overlaps = 100)
```

---

## 二、跨 Trait 图（多个 trait 做对比）

### 8. Cross-Trait 散点图 (cross_trait_scatter)

**数据来源：**

| 输入文件 | 必需字段 |
|----------|----------|
| `{trait_x}.per_gene_estimates.tsv` | `ensg`, `post_mean` |
| `{trait_y}.per_gene_estimates.tsv` | `ensg`, `post_mean` |
| `gencode_v41_gname_gid_ALL_sorted_onlyID` | 2列无表头 |

**数据处理：**
```r
# merge 两个 posterior by ensg，列名改为 effect_x / effect_y
df$gene <- label_from_ensg(df$ensg, gene_map_path)
df$label_score <- pmax(abs(effect_x), abs(effect_y))     # 取 max 排名
pca_line <- prcomp(cbind(effect_x, effect_y))$rotation    # PCA 方向线
```

**核心画图代码：**
```r
library(ggplot2); library(ggrepel)

g <- ggplot(df, aes(x = effect_x, y = effect_y))
g <- g + theme_classic(base_size = 18, base_family = "Helvetica")
g <- g + geom_point(alpha = 0.6, size = 2.5, color = "#4C78A8")
g <- g + geom_hline(yintercept = 0, linetype = "dashed", color = "grey60")
g <- g + geom_vline(xintercept = 0, linetype = "dashed", color = "grey60")
g <- g + geom_abline(slope = pca_slope, intercept = pca_intercept, color = "#B40426")  # PCA 线
g <- g + geom_text_repel(data = label_df, aes(label = label),
    size = 4.5, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.2)
g <- g + xlab(sprintf("%s posterior effect", x_label))
g <- g + ylab(sprintf("%s posterior effect", y_label))
```

---

### 9. Gene-Level QQ 图 (gene_level_qq)

**数据来源：**

| 输入文件 | 必需字段 |
|----------|----------|
| trait_targets TSV (meta 文件) | `trait_id`, `trait_stem`, `correlation_path` |
| 每个 trait 的 `*_geneRegulation_correlation.txt` | `beta_withShet`, `P_withShet` |

**数据处理：**
```r
# 对每个 trait：
df$signed_log10_p <- sign(beta_withShet) * (-log10(pmax(P_withShet, .Machine$double.xmin)))
df$expected <- signed_expected(nrow(df))  # 理论分位数
df$observed <- signed_log10_p             # 截断到 ±y_limit
```

**核心画图代码：**
```r
library(ggplot2)

g <- ggplot(summary_df, aes(x = expected, y = observed, color = trait_id))
g <- g + theme_classic(base_size = 18, base_family = "Helvetica")
g <- g + geom_point(alpha = 0.6, size = 2.2)
g <- g + geom_abline(intercept = 0, slope = 1, linetype = "dashed")
g <- g + xlab("Expected signed -log10(P)")
g <- g + ylab("Observed signed -log10(P)")
```

---

### 10. Program Heatmap (program_heatmap)

**数据来源：** 每个 trait 的 `programs_enrichment_K{k}_{trait_file}` + `regulators_enrichment_K{k}_{trait_file}`

**数据处理：**
```r
# 对每个 trait 调用 read_program_regulator_summary()
# 合并所有 trait 的 program_score / regulator_score → long 格式
# 聚类：dist + hclust(method = "complete") 对 programs 和 traits 分别聚类
ordered_traits <- colnames(mat)[hc_traits$order]
ordered_programs <- rownames(mat)[hc_programs$order]
```

**核心画图代码：**
```r
library(ggplot2)

g <- ggplot(df, aes(x = trait_id, y = Program, fill = value))
g <- g + theme_minimal(base_size = 16, base_family = "Helvetica")
g <- g + geom_tile(color = "white", linewidth = 0.2)
g <- g + scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0)
g <- g + xlab("Trait ID")
g <- g + ylab("Program")
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
```

---

### 11. Cross-Trait 相关性热图 (cross_trait_heatmap)

**数据来源：** 多个 trait 的 `{trait}.per_gene_estimates.tsv`（都用 `ensg`, `post_mean`）

**数据处理：**
```r
# 合并所有 trait 的 post_mean by ensg → 相关矩阵
cor_mat[i,j] <- cor(merged_df[[trait_i]], merged_df[[trait_j]], method = "pearson")
# 聚类排序
hc <- hclust(as.dist(1 - cor_mat), method = "complete")
```

**核心画图代码：**
```r
library(ggplot2)

g <- ggplot(plot_df, aes(x = trait_x, y = trait_y, fill = correlation))
g <- g + theme_minimal(base_size = 18, base_family = "Helvetica")
g <- g + geom_tile(color = "white", linewidth = 0.2)
g <- g + geom_text(aes(label = sprintf("%.2f", correlation)), size = 4)
g <- g + scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = c(-1, 1))
g <- g + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank())
g <- g + labs(fill = sprintf("%s r", method))
```

---

## 三、不依赖 Trait 的图（cNMF program 自身属性）

### 12. Program Top 基因条形图 (cnmf_program_top_genes)

**数据来源：**

| 输入文件 | 必需字段 |
|----------|----------|
| `cNMF_all.gene_spectra_score.k_{k}.dt_{threshold}.txt` | Programs 为行，ENSG 为列的矩阵 |
| `gencode_v41_gname_gid_ALL_sorted_onlyID` | 2列无表头 |

**核心画图代码：**
```r
library(ggplot2); library(dplyr)

# 对每个 program：取 top_n 基因按 weight 降序排列
g <- ggplot(top_genes_df, aes(x = weight, y = GENE_factor))
g <- g + theme_bw(base_size = 14, base_family = "Helvetica")
g <- g + geom_col(fill = "#3B4CC0", width = 0.7)
g <- g + facet_wrap(~ Program, scales = "free_y", ncol = min(5, n_programs))
g <- g + scale_y_discrete(labels = function(x) sub("^P[0-9]+_", "", x))
g <- g + xlab("Gene Spectra Score (Weight)")
g <- g + ylab("Gene")
```

---

### 13. Program 基因集富集气泡图 (cnmf_program_enrichment)

**数据来源：**

| 输入文件 | 必需字段 |
|----------|----------|
| `cNMF_all.gene_spectra_score.k_{k}.dt_{threshold}.txt` | Programs 为行，ENSG 为列 |
| `gencode_v41_gname_gid_ALL_sorted_onlyID` | 2列无表头 |
| `geneset_dir/{geneset}.txt` | 每行一个 gene symbol |

**数据处理：**
```r
# 每个 program × 每个 geneset：
# phyper(n_overlap - 1, n_gs, total_bg - n_gs, top_n, lower.tail = FALSE)
# FDR <- p.adjust(P, method = "BH")
```

**核心画图代码：**
```r
library(ggplot2)

g <- ggplot(enrichment_results, aes(x = Program, y = Geneset))
g <- g + theme_minimal(base_size = 14, base_family = "Helvetica")
g <- g + geom_point(aes(size = OverlapCount, color = NegLog10P))
g <- g + scale_color_gradient(low = "grey80", high = "#B40426", name = "-log10(P)")
g <- g + scale_size_continuous(name = "Overlap Genes")
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g <- g + labs(x = "Program", y = "Gene Set")
```

---

### 14. Program 间共调控散点图 (cnmf_corregulation_scatter)

**数据来源：**

| 输入文件 | 必需字段 |
|----------|----------|
| `K{k}_program{i}_perturb_effects.txt` (i=1..k) | `GENE`, `lm_es`(第2列), `lm_p`(第3列) |

**核心画图代码：**
```r
library(ggplot2)

# 筛选两个 program 中 FDR < 0.05 的基因，x=beta_a, y=beta_b
g <- ggplot(df, aes(x = x, y = y))
g <- g + theme_classic(base_size = 24, base_family = "Helvetica")
g <- g + geom_point(shape = 21, color = "black", fill = NA, stroke = 1, size = 2.5, alpha = 0.8)
g <- g + geom_vline(xintercept = 0, linetype = "dashed")
g <- g + geom_hline(yintercept = 0, linetype = "dashed")
g <- g + geom_smooth(method = "loess", method.args = list(degree = 1), span = 0.25, se = FALSE)
g <- g + xlab(paste0("Effect size on ", program_a))
g <- g + ylab(paste0("Effect size on ", program_b))
```

---

## 四、Legacy Plot 补充（`plot-*` 系列专有图）

以下图只有 legacy `plot-*` workflow 产出，`figures-*` 中没有对应实现：

### 15. GWAS vs LoF 通路富集散点图 (Fig2C)

| 输入文件 | 说明 |
|----------|------|
| `enrichment_result/{trait}_GWAS_GO.txt` | GWAS GO 富集结果 |
| `enrichment_result/{trait}_LOF_GB_GO.txt` | LoF GeneBayes GO 富集结果 |

散点图：GWAS -log10(P) vs LoF -log10(P)，标注选定通路。

### 16. 多变量回归森林图 (Fig5F)

| 输入文件 | 说明 |
|----------|------|
| cNMF regulation 矩阵 | 所有 program 的 perturb effects |
| 多个 trait 的 posterior | MCH, RDW, IRF |

水平误差棒图：标准化回归系数 + 95% CI。

### 17. 跨细胞类型组合热图 (Multi-Cell Heatmap)

| 输入文件 | 说明 |
|----------|------|
| `multiCell_model/Regulators_correlation.txt` | program 间 regulator 相关性 |
| `multiCell_model/Program_genes_Jaccard.txt` | program 间基因 Jaccard 指数 |
| 多个 trait × 多细胞类型的 enrichment | regulator enrichment 文件 |

ComplexHeatmap 双面板：上三角 = regulator 相关性 + trait 条形图，下三角 = Jaccard 指数。

### 18. 条件回归热图 (IGF-1 / HbA1c)

| 输入文件 | 说明 |
|----------|------|
| 多细胞类型 cNMF regulation | 所有程序的 perturb effects |
| 目标 trait posterior (IGF-1/HbA1c) | `*.per_gene_estimates.tsv` |

逐次回归掉 top1/2/3 regulator 后的 -log10(P) 变化的组合热图。

### 19. 交叉验证 R² 条形图 (FigS7)

| 输入文件 | 说明 |
|----------|------|
| cNMF spectra + regulation | program 基因和效应 |
| MCH/RDW/IRF posterior | 目标 trait |

80/20 交叉验证的 R² 中位数条形图，对比 3-program 选择模型 vs 随机模型 vs 单 program。

### 20. 置换检验结果图

| 输入文件 | 说明 |
|----------|------|
| `permutation_test/P{N}/{trait}_program{N}_regulator{N}_LOF{thresh}.txt` | 置换结果 |

直方图（观察值 vs 置换 null 分布）+ concordant/discordant 散点图。

---

## 共享工具函数 (helpers.R)

所有 figures R 脚本共用 `src/paper_pipeline/figures/scripts/helpers.R`，关键函数：

| 函数 | 作用 |
|------|------|
| `read_program_regulator_summary(dir, trait_file, k, label_programs)` | 读取+合并 program/regulator enrichment，计算 score 和 color |
| `read_regulation_matrices(dir, k)` | 读取所有 program 的 perturb effects，组装 beta 和 P 矩阵 |
| `parse_variant_table(gwas)` | 解析 GWAS variant (CHROM:POS)，过滤低置信度变异 |
| `add_cumulative_coordinates(gwas_df)` | 计算 GWAS 全基因组累加坐标和染色体轴刻度 |
| `read_gene_annotation_table(path)` | 读取基因注释（兼容 GTF 和简单表格） |
| `read_geneset_members(dir, geneset_name)` | 读取基因集成员列表 |
| `read_gene_map(path)` | 读取 ENSG→gene symbol 映射（2列无表头） |
| `label_from_ensg(ensg, map_path)` | 批量 ENSG→gene symbol 转换 |

所有脚本参数通过命令行 `commandArgs(trailingOnly = TRUE)` 传入，输出包括 `.tsv` 数据表和 `.pdf`/`.png` 图像。

---

## 2026-05 更新说明

### `gwas_manhattan`

- 当 GWAS 全量文件行数超过 `100000` 时，绘图前会对背景点下采样。
- 下采样目标为：
  - `50000 + ceil(1% * total_rows)`
  - 最终不超过 `300000`
- 显著点和需要标注的点优先保留，其余背景点随机抽样。

### `trait_program_gene_panel`

当前 `trait_program_gene_panel` 已按论文 Extended Data Fig. 7 的筛选语义更新。

核心参数：

- `program_n`
- `regulator_n`
- `loading_top_n`
- `random_iterations`
- `shet_path`

当前选择规则：

- 左侧 `program burden`
  - 使用 `loading_top_n` 个 loading genes（默认 `200`）
  - 对每个 program 基于 `shet` 分箱匹配的随机基因集重算 enrichment P
  - 按该 P 排序，取前 `program_n`
- 右侧 `regulator-burden`
  - 使用所有 program 的 regulator burden 与 `shet`
  - 通过 `leaps::regsubsets` 选择前 `regulator_n` 个 program

基因展示过滤：

- trait hit 需满足 `|gamma| > 0.1`（默认）
- 左侧只展示与 top loading genes 重叠的 hits
- 右侧只展示 regulator `FDR < 0.05` 且与 trait hit 重叠的 genes

注意：

- `program_n=5`、`regulator_n=3` 是当前示例默认值，对应论文中 MCH 的设置。
- 这一步是 figure 层选择逻辑；如果上游 `ProgramLevel` 仍按旧的 `top100` 预计算，图中的 program 颜色/score 可能与左侧重新筛选使用的 `top200` 依据不完全同源。
