args <- commandArgs(trailingOnly = TRUE)

posterior_path <- as.character(args[1])
limma_path      <- as.character(args[2])
shet_path       <- as.character(args[3])
trait_label     <- as.character(args[4])
table_path      <- as.character(args[5])
plot_prefix     <- as.character(args[6])
top_n_labels    <- if (length(args) >= 7) as.numeric(args[7]) else 8
highlight_genes_raw <- if (length(args) >= 8) as.character(args[8]) else ""
y_limit         <- if (length(args) >= 9) as.numeric(args[9]) else 8
gene_map_path   <- if (length(args) >= 10) as.character(args[10]) else ""

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

highlight_genes <- trimws(strsplit(highlight_genes_raw, ",", fixed = TRUE)[[1]])
highlight_genes <- highlight_genes[nzchar(highlight_genes)]

# ── 读输入 ──
options(scipen = 10)
posterior_df <- data.table::fread(posterior_path, data.table = FALSE)
if (!all(c("ensg", "post_mean") %in% names(posterior_df))) {
  stop(sprintf("Posterior file %s must have columns: ensg, post_mean", posterior_path))
}
posterior_df$ensg <- normalize_ensg_ids(posterior_df$ensg)
posterior_df$post_mean[is.infinite(posterior_df$post_mean)] <- NA

# limma logFC 矩阵（行=表达基因, 列=扰动基因）
limma_raw <- data.table::fread(limma_path, data.table = FALSE)
limma_mat <- as.matrix(limma_raw[, -1, drop = FALSE])
row.names(limma_mat) <- limma_raw[[1]]
target_ids <- row.names(limma_mat)
perturb_ids <- colnames(limma_mat)

# shet 协变量
shet_df <- utils::read.table(shet_path, header = TRUE, stringsAsFactors = FALSE)
shet_df$ensg <- normalize_ensg_ids(shet_df$ensg)

# ── 基因映射：ensg ↔ gene symbol ──
lookups <- build_gene_id_lookups(gene_map_path)
gene_lookup <- lookups$gene_lookup
target_ensg <- lookups$resolve_to_ensg(target_ids)
perturb_ensg <- lookups$resolve_to_ensg(perturb_ids)

# ── 对每个 target gene 算相关性 ──
shet_sub <- shet_df[shet_df$ensg %in% posterior_df$ensg, ]
n_targets <- nrow(limma_mat)

corr_list <- lapply(seq_len(n_targets), function(i) {
  pb <- limma_mat[i, ]  # perturb_beta 向量（扰动基因索引）
  pb_df <- data.frame(
    gene         = perturb_ids,
    perturb_beta = as.numeric(pb),
    stringsAsFactors = FALSE
  )
  pb_df$ensg <- perturb_ensg

  df <- merge(pb_df, posterior_df[, c("ensg", "post_mean")], by = "ensg", all = FALSE)
  df <- df[!is.na(df$ensg) & !is.na(df$post_mean) & !is.na(df$perturb_beta), ]
  if (!is.na(target_ensg[i])) {
    df <- df[df$ensg != target_ensg[i], ]  # 去掉靶基因本身
  }
  df <- merge(df, shet_sub, by = "ensg", all = FALSE)
  if (nrow(df) < 10) return(NULL)
  if (stats::sd(df$post_mean, na.rm = TRUE) == 0 || stats::sd(df$perturb_beta, na.rm = TRUE) == 0) {
    return(NULL)
  }

  df$post_mean_sc   <- scale(df$post_mean)
  df$perturb_beta_sc <- scale(df$perturb_beta)

  fit <- tryCatch(lm(post_mean_sc ~ perturb_beta_sc + shet, data = df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  sm  <- summary(fit)$coefficients
  if (!"perturb_beta_sc" %in% row.names(sm)) return(NULL)
  ct  <- tryCatch(cor.test(df$perturb_beta, df$post_mean, method = "pearson"),
                  error = function(e) list(p.value = NA, estimate = NA, conf.int = c(NA, NA)))

  data.frame(
    ensg                 = target_ensg[i],
    P_withShet           = sm["perturb_beta_sc", 4],
    beta_withShet        = sm["perturb_beta_sc", 1],
    betaSE_withShet      = sm["perturb_beta_sc", 2],
    P_pearson            = ct$p.value,
    R_pearson            = ct$estimate,
    R_pearson_CIlower    = ct$conf.int[1],
    R_pearson_CIupper    = ct$conf.int[2],
    stringsAsFactors = FALSE
  )
})

correlation_df <- do.call(rbind, corr_list[!vapply(corr_list, is.null, logical(1))])
if (is.null(correlation_df) || nrow(correlation_df) == 0) {
  unresolved_targets <- unique(target_ids[is.na(target_ensg)])
  unresolved_perturb <- unique(perturb_ids[is.na(perturb_ensg)])
  stop(sprintf(
    paste(
      "No gene-level scatter points produced valid correlation results",
      "(limma rows=%d, resolved target ENSG=%d, perturb columns=%d, resolved perturb ENSG=%d, posterior genes=%d, shet genes=%d).",
      "Example unresolved target IDs: %s.",
      "Example unresolved perturb IDs: %s."
    ),
    length(target_ids),
    sum(!is.na(target_ensg)),
    length(perturb_ids),
    sum(!is.na(perturb_ensg)),
    length(unique(posterior_df$ensg)),
    length(unique(shet_sub$ensg)),
    if (length(unresolved_targets) > 0) paste(head(unresolved_targets, 5), collapse = ", ") else "none",
    if (length(unresolved_perturb) > 0) paste(head(unresolved_perturb, 5), collapse = ", ") else "none"
  ))
}

# ── 合并 posterior + correlation ──
correlation_df <- correlation_df[!is.na(correlation_df$ensg) & nzchar(correlation_df$ensg), , drop = FALSE]
if (nrow(correlation_df) == 0) {
  stop("No scatter rows remained after ENSG resolution")
}
df <- merge(correlation_df, posterior_df[, c("ensg", "post_mean")], by = "ensg", all = FALSE)
df$gene <- gene_lookup[df$ensg]
df$gene[is.na(df$gene)] <- df$ensg[is.na(df$gene)]

df$signed_log10_p <- sign(df$beta_withShet) * (-log10(pmax(df$P_withShet, .Machine$double.xmin)))
df$fdr <- p.adjust(df$P_withShet, method = "BH")
df$abs_post_mean <- abs(df$post_mean)
df$abs_signed_log10_p <- abs(df$signed_log10_p)
df$post_mean_sign <- ifelse(df$post_mean > 0, "positive", ifelse(df$post_mean < 0, "negative", "zero"))
df$regulation_sign <- ifelse(df$beta_withShet > 0, "positive", ifelse(df$beta_withShet < 0, "negative", "zero"))
df$is_concordant <- df$post_mean_sign == df$regulation_sign &
  df$post_mean_sign != "zero" &
  df$regulation_sign != "zero"
df$is_discordant <- df$post_mean_sign != df$regulation_sign &
  df$post_mean_sign != "zero" &
  df$regulation_sign != "zero"
df$is_concordant[is.na(df$is_concordant)] <- FALSE
df$is_discordant[is.na(df$is_discordant)] <- FALSE
df$is_reg_sig <- !is.na(df$fdr) & df$fdr <= 0.05
df$is_high_effect <- !is.na(df$post_mean) &
  df$abs_post_mean >= stats::quantile(df$abs_post_mean, probs = 0.95, na.rm = TRUE)
df$is_high_effect[is.na(df$is_high_effect)] <- FALSE
df$combined_score <- df$abs_post_mean * df$abs_signed_log10_p
df$evidence_class <- "background"
df$evidence_class[df$is_reg_sig & df$is_concordant] <- "regulation_supported"
df$evidence_class[df$is_reg_sig & df$is_discordant] <- "direction_discordant"
df$evidence_class[!df$is_reg_sig & df$is_high_effect] <- "posterior_high"
df$evidence_class <- factor(
  df$evidence_class,
  levels = c("background", "posterior_high", "regulation_supported", "direction_discordant")
)
df$label <- ""
df$label_reason <- ""

if (length(highlight_genes) > 0) {
  highlight_mask <- df$gene %in% highlight_genes
  df$label[highlight_mask] <- df$gene[highlight_mask]
  df$label_reason[highlight_mask] <- "highlight_gene"
}

df$label_score <- df$combined_score
top_idx <- order(df$label_score, decreasing = TRUE)
if (length(top_idx) > 0) {
  top_idx <- head(top_idx, top_n_labels)
  needs_reason <- df$label[top_idx] == ""
  df$label[top_idx] <- df$gene[top_idx]
  df$label_reason[top_idx[needs_reason]] <- "top_combined_score"
}
df$label_reason[df$label != "" & df$label_reason == ""] <- "selected"

threshold_candidates <- df$P_withShet[!is.na(df$fdr) & df$fdr <= 0.05 & !is.na(df$P_withShet)]
threshold_y <- if (length(threshold_candidates) > 0) -log10(max(threshold_candidates)) else NA_real_

df$trait_label <- trait_label
df <- df[order(-df$label_score, -abs(df$signed_log10_p)), , drop = FALSE]

ensure_parent_dir(table_path)
ensure_parent_dir(paste0(plot_prefix, ".pdf"))
utils::write.table(df, table_path, row.names = FALSE, sep = "\t", quote = FALSE)

# ── 画图 ──
library(ggplot2)
library(ggrepel)

evidence_colors <- c(
  "background" = "grey78",
  "posterior_high" = "#7E8DA6",
  "regulation_supported" = "#B40426",
  "direction_discordant" = "#3B4CC0"
)
evidence_labels <- c(
  "background" = "Background",
  "posterior_high" = "High posterior effect",
  "regulation_supported" = "Concordant regulation FDR <= 0.05",
  "direction_discordant" = "Discordant regulation FDR <= 0.05"
)
evidence_shapes <- c(
  "background" = 16,
  "posterior_high" = 16,
  "regulation_supported" = 16,
  "direction_discordant" = 17
)
df$point_size <- ifelse(df$is_reg_sig, 2.9, ifelse(df$is_high_effect, 2.3, 1.7))
df$point_alpha <- ifelse(df$evidence_class == "background", 0.35, 0.82)
subtitle_text <- sprintf(
  "Each point is a gene; y sign follows beta_withShet. FDR is BH-adjusted across plotted genes (n=%s).",
  format(nrow(df), big.mark = ",", scientific = FALSE)
)
caption_text <- "Concordant genes have matching posterior and perturb-seq regulation signs; triangles mark significant discordant regulation."

g <- ggplot(df, aes(x = post_mean, y = signed_log10_p))
g <- g + theme_classic(base_size = 17, base_family = "Helvetica")
g <- g + geom_point(
  aes(color = evidence_class, shape = evidence_class, size = point_size, alpha = point_alpha),
  stroke = 0.2
)
g <- g + geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35)
g <- g + geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35)
if (!is.na(threshold_y)) {
  g <- g + geom_hline(yintercept = threshold_y, linetype = "dotted", color = "#B40426", linewidth = 0.45)
  g <- g + geom_hline(yintercept = -threshold_y, linetype = "dotted", color = "#B40426", linewidth = 0.45)
  g <- g + annotate(
    "text",
    x = Inf,
    y = threshold_y,
    label = "FDR 0.05",
    hjust = 1.04,
    vjust = -0.45,
    size = 3.6,
    color = "#B40426"
  )
}
label_df <- df[df$label != "", , drop = FALSE]
if (nrow(label_df) > 0) {
  g <- g + geom_text_repel(
    data = label_df,
    aes(label = label, color = evidence_class),
    size = 4.2,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.45,
    point.padding = 0.25,
    min.segment.length = 0,
    show.legend = FALSE
  )
}
g <- g + scale_color_manual(values = evidence_colors, labels = evidence_labels, drop = FALSE)
g <- g + scale_shape_manual(values = evidence_shapes, labels = evidence_labels, drop = FALSE)
g <- g + scale_size_identity()
g <- g + scale_alpha_identity()
g <- g + xlab(sprintf("%s posterior effect", trait_label))
g <- g + ylab("Perturb-seq regulation evidence, signed -log10(P)")
g <- g + labs(
  title = sprintf("Gene-level posterior vs perturb-seq evidence: %s", trait_label),
  subtitle = subtitle_text,
  caption = caption_text,
  color = "Evidence class",
  shape = "Evidence class"
)
g <- g + coord_cartesian(ylim = c(-y_limit, y_limit))
g <- g + theme(
  legend.position = "bottom",
  legend.title = element_text(face = "bold"),
  legend.box = "vertical",
  plot.title = element_text(face = "bold"),
  plot.subtitle = element_text(size = 11, color = "grey30"),
  plot.caption = element_text(size = 10, color = "grey35", hjust = 0),
  panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
  panel.grid.minor = element_blank()
)

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = 11, height = 10)
ggsave(plot = g, filename = paste0(plot_prefix, ".png"), width = 11, height = 10, dpi = 300)
