args <- commandArgs(trailingOnly = TRUE)

gwas_path <- as.character(args[1])
gene_annotation_path <- as.character(args[2])
source_id <- as.character(args[3])
locus_label <- as.character(args[4])
chrom <- as.character(args[5])
start_pos <- as.numeric(args[6])
end_pos <- as.numeric(args[7])
flank_bp <- as.numeric(args[8])
table_prefix <- as.character(args[9])
plot_prefix <- as.character(args[10])
genomewide_threshold <- if (length(args) >= 11) as.numeric(args[11]) else 5e-8
label_top_n <- if (length(args) >= 12) as.numeric(args[12]) else 6

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

gwas <- data.table::fread(gwas_path, data.table = FALSE)
gwas_df <- parse_variant_table(gwas)

region_start <- start_pos - flank_bp
region_end <- end_pos + flank_bp
region_df <- gwas_df[
  gwas_df$CHROM == sub("^chr", "", chrom) &
    gwas_df$POS >= region_start &
    gwas_df$POS <= region_end,
  ,
  drop = FALSE
]
if (nrow(region_df) == 0) {
  stop(sprintf("No GWAS variants found in region %s:%s-%s for %s", chrom, region_start, region_end, gwas_path))
}
region_df$neg_log10_p <- -log10(pmax(region_df$P, .Machine$double.xmin))
region_df$in_core_region <- region_df$POS >= start_pos & region_df$POS <= end_pos

gene_annotation <- read_gene_annotation_table(gene_annotation_path)

gene_df <- gene_annotation[
  sub("^chr", "", as.character(gene_annotation$chr)) == sub("^chr", "", chrom) &
    gene_annotation$end >= region_start &
    gene_annotation$start <= region_end,
  c("gene", "start", "end"),
  drop = FALSE
]
if (nrow(gene_df) > 0) {
  gene_df <- gene_df[order(gene_df$start), , drop = FALSE]
}

top_hits <- region_df[order(region_df$P, region_df$POS), , drop = FALSE]
top_hits <- head(top_hits, label_top_n)
if (nrow(top_hits) > 0) {
  top_hits$label <- paste0("chr", top_hits$CHROM, ":", format(top_hits$POS, scientific = FALSE))
}

summary_df <- data.frame(
  source_id = source_id,
  locus_label = locus_label,
  chrom = chrom,
  region_start = region_start,
  region_end = region_end,
  n_variants = nrow(region_df),
  n_core_region_variants = sum(region_df$in_core_region),
  n_genes = nrow(gene_df),
  n_top_hits = nrow(top_hits),
  stringsAsFactors = FALSE
)

ensure_parent_dir(paste0(table_prefix, "_variants.tsv"))
ensure_parent_dir(paste0(plot_prefix, ".pdf"))
utils::write.table(summary_df, paste0(table_prefix, "_summary.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
utils::write.table(region_df, paste0(table_prefix, "_variants.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
utils::write.table(gene_df, paste0(table_prefix, "_genes.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
utils::write.table(top_hits, paste0(table_prefix, "_top_hits.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

library(ggplot2)
library(ggrepel)

max_y <- max(region_df$neg_log10_p, na.rm = TRUE)
min_y <- 0
gene_track_df <- gene_df
if (nrow(gene_track_df) > 0) {
  gene_track_df$track_index <- (seq_len(nrow(gene_track_df)) - 1) %% 3
  gene_track_df$track_y <- -(gene_track_df$track_index + 1) * max(max_y * 0.08, 0.5)
  min_y <- min(gene_track_df$track_y, na.rm = TRUE) - max(max_y * 0.08, 0.5)
}

g <- ggplot(region_df, aes(x = POS, y = neg_log10_p))
g <- g + theme_classic(base_size = 18, base_family = "Helvetica")
g <- g + geom_point(aes(color = in_core_region), alpha = 0.7, size = 2.2)
g <- g + scale_color_manual(values = c("TRUE" = "#B40426", "FALSE" = "grey65"), guide = "none")
g <- g + geom_hline(yintercept = -log10(genomewide_threshold), linetype = "dashed", color = "#B40426")
g <- g + geom_vline(xintercept = c(start_pos, end_pos), linetype = "dotted", color = "grey50")
if (nrow(gene_track_df) > 0) {
  g <- g + geom_segment(
    data = gene_track_df,
    aes(x = start, xend = end, y = track_y, yend = track_y),
    inherit.aes = FALSE,
    linewidth = 1.1,
    color = "#2F4F4F"
  )
  g <- g + geom_text(
    data = gene_track_df,
    aes(x = (start + end) / 2, y = track_y, label = gene),
    inherit.aes = FALSE,
    vjust = 1.6,
    size = 3.6,
    color = "#2F4F4F"
  )
}
if (nrow(top_hits) > 0) {
  g <- g + geom_text_repel(
    data = top_hits,
    aes(label = label),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.2,
    min.segment.length = 0
  )
}
g <- g + xlab(sprintf("chr%s position", sub("^chr", "", chrom)))
g <- g + ylab("-log10(P)")
g <- g + ggtitle(locus_label)
g <- g + coord_cartesian(ylim = c(min_y, NA))

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = 12, height = 8)
ggsave(plot = g, filename = paste0(plot_prefix, ".png"), width = 12, height = 8, dpi = 300)
