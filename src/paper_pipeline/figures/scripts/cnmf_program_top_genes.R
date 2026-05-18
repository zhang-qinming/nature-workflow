args <- commandArgs(trailingOnly = TRUE)

spectra_path <- as.character(args[1])
gene_map_path <- as.character(args[2])
k <- as.numeric(args[3])
programs_raw <- as.character(args[4])
top_n <- as.numeric(args[5])
table_prefix <- as.character(args[6])
plot_prefix <- as.character(args[7])

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

programs <- parse_program_vector(programs_raw)
if (length(programs) == 0) {
  programs <- paste0("P", seq_len(k))
}

# Read spectra matrix: Programs are rows, ENSG are columns
GEP_raw <- data.table::fread(spectra_path, data.table = FALSE)
# Extract row names which might be in the first column or index
if (is.character(GEP_raw[,1])) {
    row.names(GEP_raw) <- GEP_raw[,1]
    GEP_raw <- GEP_raw[,-1]
}
GEP <- t(GEP_raw)
colnames(GEP) <- paste0("P", seq_len(k))

# Map ENSG to Symbols
gene_map <- read_gene_map(gene_map_path)
corresp <- gene_map[!duplicated(gene_map$ensg), ]
row.names(corresp) <- corresp$ensg

# Apply mapping
gene_symbols <- corresp[row.names(GEP), "gene"]
gene_symbols[is.na(gene_symbols)] <- row.names(GEP)[is.na(gene_symbols)]
GEP_df <- data.frame(GENE = gene_symbols, GEP, stringsAsFactors = FALSE)

ensure_parent_dir(paste0(table_prefix, "_top_genes.tsv"))
ensure_parent_dir(paste0(plot_prefix, "_top_genes.pdf"))

library(ggplot2)
library(dplyr)

top_genes_df <- data.frame()

for (prog in programs) {
  tmp <- GEP_df[, c("GENE", prog)]
  colnames(tmp) <- c("GENE", "weight")
  tmp <- tmp[order(tmp$weight, decreasing = TRUE), ]
  tmp_top <- head(tmp, top_n)
  tmp_top$Program <- prog
  tmp_top$rank <- seq_len(nrow(tmp_top))
  top_genes_df <- rbind(top_genes_df, tmp_top)
}

utils::write.table(top_genes_df, paste0(table_prefix, "_top_genes.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

# We will plot each program as a facet, or we can just plot them in a large grid.
# To ensure order within facets, we use tidytext-like reordering or an ordered factor
top_genes_df$Program <- factor(top_genes_df$Program, levels = programs)
top_genes_df <- top_genes_df[order(top_genes_df$Program, top_genes_df$weight), ]
top_genes_df$GENE_factor <- factor(paste(top_genes_df$Program, top_genes_df$GENE, sep="_"), levels = paste(top_genes_df$Program, top_genes_df$GENE, sep="_"))

g <- ggplot(top_genes_df, aes(x = weight, y = GENE_factor))
g <- g + theme_bw(base_size = 14, base_family = "Helvetica")
g <- g + geom_col(fill = "#3B4CC0", width=0.7)
g <- g + facet_wrap(~ Program, scales = "free_y", ncol = min(5, length(programs)))
g <- g + scale_y_discrete(labels = function(x) sub("^P[0-9]+_", "", x))
g <- g + xlab("Gene Spectra Score (Weight)")
g <- g + ylab("Gene")
g <- g + theme(
  strip.background = element_rect(fill = "grey90", color = NA),
  panel.grid.major.y = element_blank()
)

# calculate height based on rows and facets
n_rows <- ceiling(length(programs) / 5)
plot_height <- max(6, n_rows * (top_n * 0.2 + 1))
plot_width <- min(16, length(programs) * 3 + 2)

ggsave(plot = g, filename = paste0(plot_prefix, "_top_genes.pdf"), width = plot_width, height = plot_height)
ggsave(plot = g, filename = paste0(plot_prefix, "_top_genes.png"), width = plot_width, height = plot_height, dpi = 300)
