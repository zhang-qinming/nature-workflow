args <- commandArgs(trailingOnly = TRUE)

association_dir <- as.character(args[1])
trait_file <- as.character(args[2])
k <- as.numeric(args[3])
unused_arg4 <- if (length(args) >= 4) as.character(args[4]) else ""
table_path <- as.character(args[5])
plot_prefix <- as.character(args[6])

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

df <- read_program_regulator_summary(association_dir, trait_file, k, character())

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

g <- ggplot(df, aes(x = program_score, y = regulator_score, color = color))
g <- g + theme_classic(base_size = 24, base_family = "Helvetica")
g <- g + geom_point(size = 4)
g <- g + geom_vline(xintercept = 0, linetype = "dashed")
g <- g + geom_hline(yintercept = 0, linetype = "dashed")
g <- g + scale_color_manual(
  values = c(
    "other" = "grey60",
    "program_enriched" = "#FEA601",
    "regulator_enriched" = "#4783B5",
    "both_enriched" = "#34A853"
  )
)
g <- g + theme(legend.title = element_blank())
g <- g + xlab("Program burden effect, signed -log10(P)")
g <- g + ylab("Regulator-burden correlation, signed -log10(P)")

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = 14, height = 10)
