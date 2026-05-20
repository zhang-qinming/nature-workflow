args <- commandArgs(trailingOnly = TRUE)

posterior_path <- as.character(args[1])
gene_map_path <- as.character(args[2])
geneset_dir <- as.character(args[3])
highlight_genesets_raw <- if (length(args) >= 4) as.character(args[4]) else ""
genes_path <- as.character(args[5])
plot_prefix <- as.character(args[6])
label_fdr_threshold <- if (length(args) >= 7) as.numeric(args[7]) else 0.01
line_fdr_threshold <- if (length(args) >= 8) as.numeric(args[8]) else 0.1
data_genesets_raw <- if (length(args) >= 9 && nzchar(args[9])) as.character(args[9]) else ""
spectra_path <- if (length(args) >= 10 && nzchar(args[10])) as.character(args[10]) else ""
program_k <- if (length(args) >= 11) as.numeric(args[11]) else 60
top_n_program_genes <- if (length(args) >= 12) as.numeric(args[12]) else 100

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

highlight_genesets <- trimws(strsplit(highlight_genesets_raw, ",", fixed = TRUE)[[1]])
highlight_genesets <- highlight_genesets[nzchar(highlight_genesets)]
data_genesets <- if (nzchar(data_genesets_raw)) {
  sets <- trimws(strsplit(data_genesets_raw, ",", fixed = TRUE)[[1]])
  sets[nzchar(sets)]
} else {
  character()
}

posterior <- data.table::fread(posterior_path, data.table = FALSE)
required_columns <- c("ensg", "post_mean", "post2_mean", "lower_95", "upper_95")
missing_columns <- setdiff(required_columns, names(posterior))
if (length(missing_columns) > 0) {
  stop(sprintf("Posterior file %s is missing columns: %s", posterior_path, paste(missing_columns, collapse = ", ")))
}

posterior <- posterior[!is.na(posterior$post_mean) & !is.na(posterior$post2_mean), , drop = FALSE]
posterior$posterior_var <- pmax(posterior$post2_mean - posterior$post_mean^2, 0)
posterior$posterior_sd <- sqrt(posterior$posterior_var)
posterior$z <- ifelse(posterior$posterior_sd > 0, posterior$post_mean / posterior$posterior_sd, NA_real_)
posterior$P <- ifelse(
  !is.na(posterior$z),
  2 * pnorm(-abs(posterior$z)),
  NA_real_
)
posterior$P[is.na(posterior$P)] <- 1

df <- posterior[, c("ensg", "post_mean", "posterior_sd", "P", "lower_95", "upper_95"), drop = FALSE]
df$FDR <- p.adjust(df$P, method = "BH")
df$LABEL <- label_from_ensg(df$ensg, gene_map_path)

df$geneset_color <- "other"
for (geneset_name in rev(highlight_genesets)) {
  members <- read_geneset_members(geneset_dir, geneset_name)
  mask <- df$LABEL %in% members
  df$geneset_color[mask] <- friendly_geneset_name(geneset_name)
}

feature_annotations <- annotate_gene_features(
  genes = df$LABEL,
  geneset_dir = geneset_dir,
  geneset_names = data_genesets,
  spectra_path = spectra_path,
  gene_map_path = gene_map_path,
  k = as.integer(program_k),
  top_n_program_genes = as.integer(top_n_program_genes)
)
df$geneset <- feature_annotations$geneset
df$program <- feature_annotations$program
df$geneset[df$geneset == ""] <- "other"
df$program[df$program == ""] <- "other"

df$label <- ifelse(df$geneset_color != "other" & df$FDR <= label_fdr_threshold, df$LABEL, "")
df$neg_log10_p <- -log10(pmax(df$P, .Machine$double.xmin))

line_candidates <- df$P[df$FDR <= line_fdr_threshold]
threshold_y <- if (length(line_candidates) > 0) {
  -log10(max(line_candidates))
} else {
  NA_real_
}

genes_export <- data.frame(
  ensg = df$ensg,
  gene = df$LABEL,
  post_mean = df$post_mean,
  posterior_sd = df$posterior_sd,
  lower_95 = df$lower_95,
  upper_95 = df$upper_95,
  p = df$P,
  logp = df$neg_log10_p,
  fdr = df$FDR,
  geneset = df$geneset,
  program = df$program,
  stringsAsFactors = FALSE
)

hits_df <- df[df$FDR <= line_fdr_threshold, , drop = FALSE]
hits_export <- data.frame(
  ensg = hits_df$ensg,
  gene = hits_df$LABEL,
  post_mean = hits_df$post_mean,
  posterior_sd = hits_df$posterior_sd,
  lower_95 = hits_df$lower_95,
  upper_95 = hits_df$upper_95,
  p = hits_df$P,
  logp = hits_df$neg_log10_p,
  fdr = hits_df$FDR,
  geneset = hits_df$geneset,
  program = hits_df$program,
  stringsAsFactors = FALSE
)

hits_path <- gsub("_genes\\.tsv$", "_hits.tsv", genes_path)
if (hits_path == genes_path) {
  hits_path <- gsub("\\.tsv$", "_hits.tsv", genes_path)
}
if (hits_path == genes_path) {
  hits_path <- paste0(genes_path, "_hits.tsv")
}

ensure_parent_dir(genes_path)
ensure_parent_dir(paste0(plot_prefix, ".pdf"))
write.table(genes_export, genes_path, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(hits_export, hits_path, row.names = FALSE, sep = "\t", quote = FALSE)

library(ggplot2)
library(ggrepel)

palette_values <- c("other" = "grey70")
highlight_levels <- unique(df$geneset_color[df$geneset_color != "other"])
if (length(highlight_levels) > 0) {
  colors <- c("#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD")
  for (i in seq_along(highlight_levels)) {
    palette_values[highlight_levels[i]] <- colors[((i - 1) %% length(colors)) + 1]
  }
}

highlight_df <- df[df$geneset_color != "other" & df$FDR <= line_fdr_threshold, , drop = FALSE]

g <- ggplot(df, aes(x = post_mean, y = neg_log10_p, color = geneset_color))
g <- g + theme_classic(base_size = 20, base_family = "Helvetica")
g <- g + geom_point(alpha = 0.45, size = 2.5)
if (nrow(highlight_df) > 0) {
  g <- g + geom_point(data = highlight_df, alpha = 0.9, size = 3)
  g <- g + geom_text_repel(
    data = highlight_df[highlight_df$label != "", , drop = FALSE],
    aes(label = label),
    size = 5,
    max.overlaps = 100,
    show.legend = FALSE
  )
}
if (!is.na(threshold_y)) {
  g <- g + geom_hline(yintercept = threshold_y, linetype = "dashed")
}
g <- g + scale_color_manual(values = palette_values)
g <- g + theme(legend.title = element_blank())
g <- g + xlab("Posterior effect size")
g <- g + ylab("-log10(P)")

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = 14, height = 12, dpi = 300)
ggsave(plot = g, filename = paste0(plot_prefix, ".png"), width = 14, height = 12, dpi = 300)
