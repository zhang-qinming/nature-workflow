args <- commandArgs(trailingOnly = TRUE)

association_dir <- as.character(args[1])
trait_file <- as.character(args[2])
trait_id <- as.character(args[3])
k <- as.numeric(args[4])
top_n <- if (length(args) >= 5) as.numeric(args[5]) else 12
table_prefix <- as.character(args[6])
plot_prefix <- as.character(args[7])

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

df <- read_program_regulator_summary(association_dir, trait_file, k, character())
df$trait_id <- trait_id
df$program_rank_abs <- rank(-abs(df$program_score), ties.method = "first")
df$regulator_rank_abs <- rank(-abs(df$regulator_score), ties.method = "first")
df <- df[order(df$program_rank_abs, df$regulator_rank_abs), , drop = FALSE]

ensure_parent_dir(paste0(table_prefix, "_summary.tsv"))
ensure_parent_dir(paste0(plot_prefix, "_program_score.pdf"))
utils::write.table(df, paste0(table_prefix, "_summary.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

library(ggplot2)

plot_top_metric <- function(metric_name, rank_name, title_text, file_stem) {
  tmp <- df[order(df[[rank_name]]), c("Program", metric_name, rank_name), drop = FALSE]
  tmp <- head(tmp, top_n)
  tmp$Program <- factor(tmp$Program, levels = rev(tmp$Program))
  colnames(tmp)[2] <- "value"

  g <- ggplot(tmp, aes(x = Program, y = value, fill = value > 0))
  g <- g + theme_classic(base_size = 18, base_family = "Helvetica")
  g <- g + geom_col(width = 0.75)
  g <- g + coord_flip()
  g <- g + geom_hline(yintercept = 0, linetype = "dashed")
  g <- g + scale_fill_manual(values = c("TRUE" = "#B40426", "FALSE" = "#3B4CC0"), guide = "none")
  g <- g + xlab("Program")
  g <- g + ylab(title_text)

  ggsave(plot = g, filename = paste0(plot_prefix, "_", file_stem, ".pdf"), width = 10, height = 8)
  ggsave(plot = g, filename = paste0(plot_prefix, "_", file_stem, ".png"), width = 10, height = 8, dpi = 300)
}

plot_top_metric("program_score", "program_rank_abs", "Program burden effect, signed -log10(P)", "program_score")
plot_top_metric("regulator_score", "regulator_rank_abs", "Regulator-burden correlation, signed -log10(P)", "regulator_score")
