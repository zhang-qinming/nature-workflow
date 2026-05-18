args <- commandArgs(trailingOnly = TRUE)

posterior_path_x <- as.character(args[1])
posterior_path_y <- as.character(args[2])
x_label <- as.character(args[3])
y_label <- as.character(args[4])
gene_map_path <- as.character(args[5])
table_path <- as.character(args[6])
plot_prefix <- as.character(args[7])
top_n_labels <- if (length(args) >= 8) as.numeric(args[8]) else 12
highlight_genes_raw <- if (length(args) >= 9) as.character(args[9]) else ""

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

highlight_genes <- trimws(strsplit(highlight_genes_raw, ",", fixed = TRUE)[[1]])
highlight_genes <- highlight_genes[nzchar(highlight_genes)]

x_df <- data.table::fread(posterior_path_x, data.table = FALSE)
y_df <- data.table::fread(posterior_path_y, data.table = FALSE)

required_columns <- c("ensg", "post_mean")
missing_x <- setdiff(required_columns, names(x_df))
missing_y <- setdiff(required_columns, names(y_df))
if (length(missing_x) > 0) {
  stop(sprintf("Posterior file %s is missing columns: %s", posterior_path_x, paste(missing_x, collapse = ", ")))
}
if (length(missing_y) > 0) {
  stop(sprintf("Posterior file %s is missing columns: %s", posterior_path_y, paste(missing_y, collapse = ", ")))
}

colnames(x_df)[colnames(x_df) == "post_mean"] <- "effect_x"
colnames(y_df)[colnames(y_df) == "post_mean"] <- "effect_y"
df <- merge(
  x_df[, c("ensg", "effect_x")],
  y_df[, c("ensg", "effect_y")],
  by = "ensg",
  all = FALSE
)
df$gene <- label_from_ensg(df$ensg, gene_map_path)
df$label <- ""

if (length(highlight_genes) > 0) {
  df$label[df$gene %in% highlight_genes] <- df$gene[df$gene %in% highlight_genes]
}

df$label_score <- pmax(abs(df$effect_x), abs(df$effect_y))
top_idx <- order(df$label_score, decreasing = TRUE)
if (length(top_idx) > 0) {
  top_idx <- head(top_idx, top_n_labels)
  df$label[top_idx] <- df$gene[top_idx]
}

df <- df[order(-df$label_score), , drop = FALSE]

ensure_parent_dir(table_path)
ensure_parent_dir(paste0(plot_prefix, ".pdf"))
utils::write.table(df, table_path, row.names = FALSE, sep = "\t", quote = FALSE)

library(ggplot2)
library(ggrepel)

pca_line <- tryCatch({
  rotation <- stats::prcomp(cbind(df$effect_x, df$effect_y))$rotation
  slope <- rotation[2, 1] / rotation[1, 1]
  intercept <- mean(df$effect_y) - slope * mean(df$effect_x)
  list(slope = slope, intercept = intercept)
}, error = function(...) {
  NULL
})

g <- ggplot(df, aes(x = effect_x, y = effect_y))
g <- g + theme_classic(base_size = 18, base_family = "Helvetica")
g <- g + geom_point(alpha = 0.6, size = 2.5, color = "#4C78A8")
g <- g + geom_hline(yintercept = 0, linetype = "dashed", color = "grey60")
g <- g + geom_vline(xintercept = 0, linetype = "dashed", color = "grey60")
if (!is.null(pca_line)) {
  g <- g + geom_abline(slope = pca_line$slope, intercept = pca_line$intercept, color = "#B40426")
}
label_df <- df[df$label != "", , drop = FALSE]
if (nrow(label_df) > 0) {
  g <- g + geom_text_repel(
    data = label_df,
    aes(label = label),
    size = 4.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    min.segment.length = 0
  )
}
g <- g + xlab(sprintf("%s posterior effect", x_label))
g <- g + ylab(sprintf("%s posterior effect", y_label))

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = 11, height = 10)
ggsave(plot = g, filename = paste0(plot_prefix, ".png"), width = 11, height = 10, dpi = 300)
