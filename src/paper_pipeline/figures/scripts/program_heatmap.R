args <- commandArgs(trailingOnly = TRUE)

association_dir <- as.character(args[1])
k <- as.numeric(args[2])
trait_target_path <- as.character(args[3])
table_prefix <- as.character(args[4])
plot_prefix <- as.character(args[5])
metrics_raw <- if (length(args) >= 6) as.character(args[6]) else "program_score,regulator_score"

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

metrics <- trimws(strsplit(metrics_raw, ",", fixed = TRUE)[[1]])
metrics <- metrics[nzchar(metrics)]
if (length(metrics) == 0) {
  stop("program_heatmap.R requires at least one metric")
}

targets <- utils::read.table(trait_target_path, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
required_columns <- c("trait_id", "trait_file")
missing_columns <- setdiff(required_columns, names(targets))
if (length(missing_columns) > 0) {
  stop(sprintf("Trait target file %s is missing columns: %s", trait_target_path, paste(missing_columns, collapse = ", ")))
}

summary_long <- data.frame()
for (i in seq_len(nrow(targets))) {
  target_id <- as.character(targets$trait_id[i])
  trait_file <- as.character(targets$trait_file[i])
  df <- read_program_regulator_summary(association_dir, trait_file, k, character())
  df$trait_id <- target_id
  df$trait_file <- trait_file
  summary_long <- rbind(summary_long, df[, c("Program", "trait_id", "trait_file", "program_score", "regulator_score")])
}

summary_long$Program <- factor(summary_long$Program, levels = paste0("P", seq_len(k)))
summary_long$trait_id <- factor(summary_long$trait_id, levels = targets$trait_id)

ensure_parent_dir(paste0(table_prefix, "_long.tsv"))
ensure_parent_dir(paste0(plot_prefix, "_program_score.pdf"))
utils::write.table(summary_long, paste0(table_prefix, "_long.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

write_metric_matrix <- function(metric_name) {
  tmp <- summary_long[, c("Program", "trait_id", metric_name)]
  mat <- reshape(
    tmp,
    idvar = "Program",
    timevar = "trait_id",
    direction = "wide"
  )
  utils::write.table(mat, paste0(table_prefix, "_", metric_name, "_matrix.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
}

library(ggplot2)

plot_metric <- function(metric_name, title_text) {
  df <- summary_long[, c("Program", "trait_id", metric_name)]
  colnames(df)[3] <- "value"
  
  # Reshape for clustering
  mat <- reshape(
    df,
    idvar = "Program",
    timevar = "trait_id",
    direction = "wide"
  )
  rownames(mat) <- mat$Program
  mat$Program <- NULL
  
  # Fill NAs with 0 for clustering purposes
  mat[is.na(mat)] <- 0
  
  # Cluster traits (columns)
  dist_traits <- stats::dist(t(mat))
  hc_traits <- stats::hclust(dist_traits, method = "complete")
  ordered_traits <- colnames(mat)[hc_traits$order]
  ordered_traits <- sub("^value\\.", "", ordered_traits)
  
  # Cluster programs (rows)
  dist_programs <- stats::dist(mat)
  hc_programs <- stats::hclust(dist_programs, method = "complete")
  ordered_programs <- rownames(mat)[hc_programs$order]

  df$trait_id <- factor(df$trait_id, levels = ordered_traits)
  df$Program <- factor(df$Program, levels = ordered_programs)

  g <- ggplot(df, aes(x = trait_id, y = Program, fill = value))
  g <- g + theme_minimal(base_size = 16, base_family = "Helvetica")
  g <- g + geom_tile(color = "white", linewidth = 0.2)
  g <- g + scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0)
  g <- g + xlab("Trait ID")
  g <- g + ylab("Program")
  g <- g + labs(fill = title_text)
  g <- g + theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
  ggsave(plot = g, filename = paste0(plot_prefix, "_", metric_name, ".pdf"), width = max(8, 1.2 * nrow(targets)), height = 14)
  ggsave(plot = g, filename = paste0(plot_prefix, "_", metric_name, ".png"), width = max(8, 1.2 * nrow(targets)), height = 14, dpi = 300)
}

for (metric_name in metrics) {
  if (!(metric_name %in% c("program_score", "regulator_score"))) {
    stop(sprintf("Unsupported metric '%s'", metric_name))
  }
  write_metric_matrix(metric_name)
  title_text <- if (metric_name == "program_score") {
    "Program score"
  } else {
    "Regulator score"
  }
  plot_metric(metric_name, title_text)
}
