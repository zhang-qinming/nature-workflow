args <- commandArgs(trailingOnly = TRUE)

regulation_dir <- as.character(args[1])
k <- as.numeric(args[2])
program_a <- as.character(args[3])
program_b <- as.character(args[4])
table_path <- as.character(args[5])
plot_prefix <- as.character(args[6])

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

matrices <- read_regulation_matrices(regulation_dir, k)
beta_df <- matrices$beta
p_df <- matrices$p

target <- unique(c(
  p_df$GENE[p.adjust(p_df[[program_a]], method = "BH") < 0.05],
  p_df$GENE[p.adjust(p_df[[program_b]], method = "BH") < 0.05]
))
df <- beta_df[is.element(beta_df$GENE, target), c("GENE", program_a, program_b)]
colnames(df) <- c("GENE", "x", "y")

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

g <- ggplot(df, aes(x = x, y = y))
g <- g + theme_classic(base_size = 24, base_family = "Helvetica")
g <- g + geom_point(shape = 21, color = "black", fill = NA, stroke = 1, size = 2.5, alpha = 0.8)
g <- g + geom_vline(xintercept = 0, linetype = "dashed")
g <- g + geom_hline(yintercept = 0, linetype = "dashed")
g <- g + geom_smooth(method = "loess", method.args = list(degree = 1), span = 0.25, se = FALSE)
g <- g + theme(legend.position = "none")
g <- g + xlab(paste0("Effect size on ", program_a))
g <- g + ylab(paste0("Effect size on ", program_b))

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = 8, height = 8)
