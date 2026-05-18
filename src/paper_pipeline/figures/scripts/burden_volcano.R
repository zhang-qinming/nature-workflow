args <- commandArgs(trailingOnly = TRUE)

burden_path <- as.character(args[1])
gene_map_path <- as.character(args[2])
geneset_dir <- as.character(args[3])
highlight_genesets_raw <- if (length(args) >= 4) as.character(args[4]) else ""
table_path <- as.character(args[5])
plot_prefix <- as.character(args[6])
label_fdr_threshold <- if (length(args) >= 7) as.numeric(args[7]) else 0.01
line_fdr_threshold <- if (length(args) >= 8) as.numeric(args[8]) else 0.1

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

highlight_genesets <- trimws(strsplit(highlight_genesets_raw, ",", fixed = TRUE)[[1]])
highlight_genesets <- highlight_genesets[nzchar(highlight_genesets)]

lof <- data.table::fread(burden_path, data.table = FALSE)
required_columns <- c("beta", "standard_error", "ensg")
missing_columns <- setdiff(required_columns, names(lof))
if (length(missing_columns) > 0) {
  stop(sprintf("Burden file %s is missing columns: %s", burden_path, paste(missing_columns, collapse = ", ")))
}

lof <- data.frame(
  lof,
  P = 2 * pnorm(-abs(lof$beta / lof$standard_error)),
  stringsAsFactors = FALSE
)

df <- lof[, c("beta", "P", "ensg")]
df$FDR <- p.adjust(df$P, method = "BH")
df$LABEL <- label_from_ensg(df$ensg, gene_map_path)

# 累积全部匹配的基因集（; 分隔），同时记录第一个匹配的用于画图着色
df$geneset       <- ""
df$geneset_color <- "other"
for (geneset_name in rev(highlight_genesets)) {
  members <- read_geneset_members(geneset_dir, geneset_name)
  mask <- df$LABEL %in% members
  display <- friendly_geneset_name(geneset_name)
  # 累积全部匹配
  df$geneset[mask] <- ifelse(
    nchar(df$geneset[mask]) == 0,
    display,
    paste(df$geneset[mask], display, sep = ";")
  )
  # 画图着色用第一个（循环最后一次生效的）
  df$geneset_color[mask] <- display
}
df$geneset[df$geneset == ""] <- "other"

df$label <- ifelse(df$geneset_color != "other" & df$FDR <= label_fdr_threshold, df$LABEL, "")

df$neg_log10_p <- -log10(pmax(df$P, .Machine$double.xmin))
line_candidates <- df$P[df$FDR <= line_fdr_threshold]
threshold_y <- if (length(line_candidates) > 0) {
  -log10(max(line_candidates))
} else {
  NA_real_
}

ensure_parent_dir(table_path)
ensure_parent_dir(paste0(plot_prefix, ".pdf"))
write.table(
  df,
  table_path,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

library(ggplot2)
library(ggrepel)

palette_values <- c("other" = "grey70")
for (i in seq_along(unique(df$geneset_color[df$geneset_color != "other"]))) {
  palette_values[unique(df$geneset_color[df$geneset_color != "other"])[i]] <- c("#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD")[((i - 1) %% 5) + 1]
}

highlight_df <- df[df$geneset_color != "other" & df$FDR <= line_fdr_threshold, , drop = FALSE]

g <- ggplot(df, aes(x = beta, y = neg_log10_p, color = geneset_color))
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
g <- g + xlab("Effect size")
g <- g + ylab("-log10(P)")

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = 14, height = 12, dpi = 300)
ggsave(plot = g, filename = paste0(plot_prefix, ".png"), width = 14, height = 12, dpi = 300)
