args <- commandArgs(trailingOnly = TRUE)

target_path <- as.character(args[1])
gene_map_path <- as.character(args[2])
table_prefix <- as.character(args[3])
plot_prefix <- as.character(args[4])
method <- if (length(args) >= 5) as.character(args[5]) else "pearson"

targets <- utils::read.table(target_path, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
required_columns <- c("trait_id", "posterior_path")
missing_columns <- setdiff(required_columns, names(targets))
if (length(missing_columns) > 0) {
  stop(sprintf("cross_trait_heatmap.R target file is missing columns: %s", paste(missing_columns, collapse = ", ")))
}

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

posterior_list <- list()
for (i in seq_len(nrow(targets))) {
  trait_id <- as.character(targets$trait_id[i])
  posterior_path <- as.character(targets$posterior_path[i])
  df <- data.table::fread(posterior_path, data.table = FALSE)
  required_post <- c("ensg", "post_mean")
  missing_post <- setdiff(required_post, names(df))
  if (length(missing_post) > 0) {
    stop(sprintf("Posterior file %s is missing columns: %s", posterior_path, paste(missing_post, collapse = ", ")))
  }
  colnames(df)[colnames(df) == "post_mean"] <- trait_id
  posterior_list[[trait_id]] <- df[, c("ensg", trait_id)]
}

merged_df <- posterior_list[[1]]
if (length(posterior_list) > 1) {
  for (i in 2:length(posterior_list)) {
    merged_df <- merge(merged_df, posterior_list[[i]], by = "ensg", all = FALSE)
  }
}

trait_ids <- as.character(targets$trait_id)
cor_mat <- matrix(NA_real_, nrow = length(trait_ids), ncol = length(trait_ids), dimnames = list(trait_ids, trait_ids))
pairwise_df <- data.frame()

for (i in seq_along(trait_ids)) {
  for (j in seq_along(trait_ids)) {
    x <- merged_df[[trait_ids[i]]]
    y <- merged_df[[trait_ids[j]]]
    cor_val <- suppressWarnings(stats::cor(x, y, method = method, use = "pairwise.complete.obs"))
    cor_mat[i, j] <- cor_val
    if (j >= i) {
      pairwise_df <- rbind(
        pairwise_df,
        data.frame(
          trait_x = trait_ids[i],
          trait_y = trait_ids[j],
          correlation = cor_val,
          n_genes = sum(stats::complete.cases(x, y)),
          method = method,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

matrix_df <- data.frame(trait_id = rownames(cor_mat), cor_mat, check.names = FALSE, stringsAsFactors = FALSE)

ensure_parent_dir(paste0(table_prefix, "_matrix.tsv"))
ensure_parent_dir(paste0(plot_prefix, ".pdf"))
utils::write.table(pairwise_df, paste0(table_prefix, "_pairs.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
utils::write.table(matrix_df, paste0(table_prefix, "_matrix.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

plot_df <- data.frame()
for (i in seq_along(trait_ids)) {
  for (j in seq_along(trait_ids)) {
    plot_df <- rbind(
      plot_df,
      data.frame(
        trait_x = trait_ids[i],
        trait_y = trait_ids[j],
        correlation = cor_mat[i, j],
        stringsAsFactors = FALSE
      )
    )
  }
}
dist_mat <- stats::as.dist(1 - cor_mat)
hc <- stats::hclust(dist_mat, method = "complete")
ordered_traits <- trait_ids[hc$order]

plot_df$trait_x <- factor(plot_df$trait_x, levels = ordered_traits)
plot_df$trait_y <- factor(plot_df$trait_y, levels = rev(ordered_traits))

library(ggplot2)

g <- ggplot(plot_df, aes(x = trait_x, y = trait_y, fill = correlation))
g <- g + theme_minimal(base_size = 18, base_family = "Helvetica")
g <- g + geom_tile(color = "white", linewidth = 0.2)
g <- g + geom_text(aes(label = sprintf("%.2f", correlation)), size = 4)
g <- g + scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = c(-1, 1))
g <- g + theme(
  axis.title = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1),
  panel.grid = element_blank()
)
g <- g + labs(fill = sprintf("%s r", method))

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = max(8, 1.1 * length(trait_ids)), height = max(8, 1.0 * length(trait_ids)))
ggsave(plot = g, filename = paste0(plot_prefix, ".png"), width = max(8, 1.1 * length(trait_ids)), height = max(8, 1.0 * length(trait_ids)), dpi = 300)
