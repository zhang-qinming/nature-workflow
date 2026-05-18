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
posterior_df$post_mean[is.infinite(posterior_df$post_mean)] <- NA

# limma logFC 矩阵（行=表达基因, 列=扰动基因）
limma_raw <- data.table::fread(limma_path, data.table = FALSE)
limma_mat <- as.matrix(limma_raw[, -1, drop = FALSE])
row.names(limma_mat) <- limma_raw[[1]]
perturb_genes <- colnames(limma_mat)

# shet 协变量
shet_df <- utils::read.table(shet_path, header = TRUE, stringsAsFactors = FALSE)

# ── 基因映射：ensg ↔ gene symbol ──
gene_map <- read_gene_map(gene_map_path)
gene_lookup <- stats::setNames(gene_map$gene, gene_map$ensg)  # ensg → symbol
ensg_lookup <- stats::setNames(gene_map$ensg, gene_map$gene)  # symbol → ensg

# ── 对每个扰动基因算相关性 ──
shet_sub <- shet_df[shet_df$ensg %in% posterior_df$ensg, ]
n_genes  <- length(perturb_genes)

corr_list <- lapply(seq_len(n_genes), function(i) {
  target_gene <- perturb_genes[i]
  pb <- limma_mat[, i]  # perturb_beta 向量（基因名索引）
  pb_df <- data.frame(
    gene         = names(pb),
    perturb_beta = as.numeric(pb),
    stringsAsFactors = FALSE
  )
  pb_df$ensg <- ensg_lookup[pb_df$gene]

  df <- merge(pb_df, posterior_df[, c("ensg", "post_mean")], by = "ensg", all = FALSE)
  df <- df[!is.na(df$post_mean) & !is.na(df$perturb_beta), ]
  df <- df[df$ensg != ensg_lookup[target_gene], ]  # 去掉靶基因本身
  df <- merge(df, shet_sub, by = "ensg", all = FALSE)
  if (nrow(df) < 10) return(NULL)

  df$post_mean_sc   <- scale(df$post_mean)
  df$perturb_beta_sc <- scale(df$perturb_beta)

  fit <- lm(post_mean_sc ~ perturb_beta_sc + shet, data = df)
  sm  <- summary(fit)$coefficients
  ct  <- tryCatch(cor.test(df$perturb_beta, df$post_mean, method = "pearson"),
                  error = function(e) list(p.value = NA, estimate = NA, conf.int = c(NA, NA)))

  data.frame(
    ensg                 = ensg_lookup[target_gene],
    P_withShet           = sm[2, 4],
    beta_withShet        = sm[2, 1],
    betaSE_withShet      = sm[2, 2],
    P_pearson            = ct$p.value,
    R_pearson            = ct$estimate,
    R_pearson_CIlower    = ct$conf.int[1],
    R_pearson_CIupper    = ct$conf.int[2],
    stringsAsFactors = FALSE
  )
})

correlation_df <- do.call(rbind, corr_list[!vapply(corr_list, is.null, logical(1))])
if (is.null(correlation_df) || nrow(correlation_df) == 0) {
  stop("No perturb genes produced valid correlation results")
}

# ── 合并 posterior + correlation ──
df <- merge(correlation_df, posterior_df[, c("ensg", "post_mean")], by = "ensg", all = FALSE)
df$gene <- gene_lookup[df$ensg]
df$gene[is.na(df$gene)] <- df$ensg[is.na(df$gene)]

df$signed_log10_p <- sign(df$beta_withShet) * (-log10(pmax(df$P_withShet, .Machine$double.xmin)))
df$fdr <- p.adjust(df$P_withShet, method = "BH")
df$label <- ""

if (length(highlight_genes) > 0) {
  df$label[df$gene %in% highlight_genes] <- df$gene[df$gene %in% highlight_genes]
}

df$label_score <- abs(df$post_mean) * abs(df$signed_log10_p)
top_idx <- order(df$label_score, decreasing = TRUE)
if (length(top_idx) > 0) {
  top_idx <- head(top_idx, top_n_labels)
  df$label[top_idx] <- df$gene[top_idx]
}

threshold_candidates <- df$P_withShet[df$fdr <= 0.05]
threshold_y <- if (length(threshold_candidates) > 0) -log10(max(threshold_candidates)) else NA_real_

df$trait_label <- trait_label
df <- df[order(-df$label_score, -abs(df$signed_log10_p)), , drop = FALSE]

ensure_parent_dir(table_path)
ensure_parent_dir(paste0(plot_prefix, ".pdf"))
utils::write.table(df, table_path, row.names = FALSE, sep = "\t", quote = FALSE)

# ── 画图 ──
library(ggplot2)
library(ggrepel)

g <- ggplot(df, aes(x = post_mean, y = signed_log10_p))
g <- g + theme_classic(base_size = 18, base_family = "Helvetica")
g <- g + geom_point(alpha = 0.6, size = 2.5, color = "#4C78A8")
g <- g + geom_hline(yintercept = 0, linetype = "dashed", color = "grey60")
g <- g + geom_vline(xintercept = 0, linetype = "dashed", color = "grey60")
if (!is.na(threshold_y)) {
  g <- g + geom_hline(yintercept = threshold_y, linetype = "dotted", color = "#B40426")
  g <- g + geom_hline(yintercept = -threshold_y, linetype = "dotted", color = "#B40426")
}
label_df <- df[df$label != "", , drop = FALSE]
if (nrow(label_df) > 0) {
  g <- g + geom_text_repel(data = label_df, aes(label = label),
    size = 4.5, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.2, min.segment.length = 0)
}
g <- g + xlab(sprintf("%s posterior effect", trait_label))
g <- g + ylab("Gene-level signed -log10(P)")
g <- g + coord_cartesian(ylim = c(-y_limit, y_limit))

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = 11, height = 10)
ggsave(plot = g, filename = paste0(plot_prefix, ".png"), width = 11, height = 10, dpi = 300)
