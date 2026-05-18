args <- commandArgs(trailingOnly = TRUE)

gwas_path <- as.character(args[1])
gene_annotation_path <- as.character(args[2])
geneset_dir <- as.character(args[3])
highlight_genesets_raw <- if (length(args) >= 4) as.character(args[4]) else ""
table_path <- as.character(args[5])
plot_prefix <- as.character(args[6])
flank_bp <- if (length(args) >= 7) as.numeric(args[7]) else 50000
label_p_threshold <- if (length(args) >= 8) as.numeric(args[8]) else 1e-30
genomewide_threshold <- if (length(args) >= 9) as.numeric(args[9]) else 5e-8

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

highlight_genesets <- trimws(strsplit(highlight_genesets_raw, ",", fixed = TRUE)[[1]])
highlight_genesets <- highlight_genesets[nzchar(highlight_genesets)]

gwas <- data.table::fread(gwas_path, data.table = FALSE)
gwas_df <- parse_variant_table(gwas)
coord_data <- add_cumulative_coordinates(gwas_df)
gwas_df <- coord_data$data
axis_df <- coord_data$axis

gene_annotation <- read_gene_annotation_table(gene_annotation_path)

# 构建基因集区间表
region_rows <- list()
for (geneset_name in highlight_genesets) {
  genes <- read_geneset_members(geneset_dir, geneset_name)
  gene_rows <- gene_annotation[gene_annotation$gene %in% genes, c("gene", "chr", "start", "end")]
  if (nrow(gene_rows) == 0) next
  gene_rows$chr <- sub("^chr", "", as.character(gene_rows$chr))
  gene_rows$region_start <- gene_rows$start - flank_bp
  gene_rows$region_end   <- gene_rows$end   + flank_bp
  gene_rows$geneset_name <- geneset_name
  gene_rows$display_name <- friendly_geneset_name(geneset_name)
  region_rows[[length(region_rows) + 1]] <- gene_rows
}
region_df <- do.call(rbind, region_rows)
has_genesets <- !is.null(region_df) && nrow(region_df) > 0

# 给每个变异点累积全部匹配的基因集（; 分隔），同时记录第一个用于着色
gwas_df$geneset       <- ""
gwas_df$geneset_color <- "background"
label_candidates <- data.frame()

if (has_genesets) {
  all_chroms <- unique(gwas_df$CHROM)

  for (chr in all_chroms) {
    chr_mask <- gwas_df$CHROM == chr
    chr_regions <- region_df[region_df$chr == chr, , drop = FALSE]
    if (nrow(chr_regions) == 0) next

    # rev() 让后遍历的基因集先出现在列表前部
    chr_regions <- chr_regions[nrow(chr_regions):1, ]

    chr_pos <- gwas_df$POS[chr_mask]
    chr_geneset       <- rep("", sum(chr_mask))
    chr_geneset_color <- rep("background", sum(chr_mask))

    for (j in seq_len(nrow(chr_regions))) {
      in_region <- chr_pos >= chr_regions$region_start[j] &
                   chr_pos <= chr_regions$region_end[j]
      # 累积全部匹配
      chr_geneset[in_region] <- ifelse(
        nchar(chr_geneset[in_region]) == 0,
        chr_regions$display_name[j],
        paste(chr_geneset[in_region], chr_regions$display_name[j], sep = ";")
      )
      # 着色用最后一个 rev 后匹配到的
      chr_geneset_color[in_region] <- chr_regions$display_name[j]
    }
    gwas_df$geneset[chr_mask]       <- chr_geneset
    gwas_df$geneset_color[chr_mask] <- chr_geneset_color

    # 收集可标注位点
    sig_chr_mask <- chr_mask & gwas_df$P <= label_p_threshold & gwas_df$geneset != ""
    if (any(sig_chr_mask)) {
      label_candidates <- rbind(
        label_candidates,
        gwas_df[sig_chr_mask, c("CHROM", "POS", "P", "COORD", "geneset", "geneset_color"), drop = FALSE]
      )
    }
  }

  # 给显著位点标注最近基因名
  if (nrow(label_candidates) > 0) {
    label_candidates$LABEL <- ""
    for (j in seq_len(nrow(label_candidates))) {
      chr_j    <- label_candidates$CHROM[j]
      pos_j    <- label_candidates$POS[j]
      color_j  <- label_candidates$geneset_color[j]

      geneset_orig <- region_df$display_name == color_j & region_df$chr == chr_j
      nearby <- region_df[geneset_orig &
                          pos_j >= region_df$region_start &
                          pos_j <= region_df$region_end, , drop = FALSE]
      if (nrow(nearby) > 0) {
        label_candidates$LABEL[j] <- nearby$gene[1]
      }
    }
    label_candidates <- label_candidates[label_candidates$LABEL != "", , drop = FALSE]
    label_candidates <- label_candidates[!duplicated(label_candidates$COORD), , drop = FALSE]
  }
} else {
  label_candidates <- data.frame(
    CHROM = character(), POS = numeric(), P = numeric(),
    COORD = numeric(), geneset = character(), geneset_color = character(),
    LABEL = character(), stringsAsFactors = FALSE
  )
}

gwas_df$geneset[gwas_df$geneset == ""] <- "background"
gwas_df$neg_log10_p <- -log10(pmax(gwas_df$P, .Machine$double.xmin))
ensure_parent_dir(table_path)
ensure_parent_dir(paste0(plot_prefix, ".pdf"))

# 保存主表（geneset = 全量 ;分隔, geneset_color = 着色用）
write.table(
  gwas_df[, c("CHROM", "POS", "P", "COORD", "neg_log10_p", "geneset", "geneset_color")],
  table_path,
  row.names = FALSE, sep = "\t", quote = FALSE
)
# 显著位点单独存
hits_path <- gsub("\\.tsv$", "_hits.tsv", table_path)
if (hits_path == table_path) hits_path <- paste0(table_path, "_hits.tsv")
write.table(label_candidates, hits_path, row.names = FALSE, sep = "\t", quote = FALSE)

# 画图（用 geneset_color 着色）
library(ggplot2)
library(ggrepel)
library(ggrastr)

unique_sets <- unique(gwas_df$geneset_color[gwas_df$geneset_color != "background"])
palette_values <- c("background" = "grey80")
if (length(unique_sets) > 0) {
  highlight_colors <- c(
    "#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"
  )
  for (i in seq_along(unique_sets)) {
    palette_values[unique_sets[i]] <- highlight_colors[((i - 1) %% length(highlight_colors)) + 1]
  }
}

g <- ggplot(gwas_df, aes(x = COORD, y = neg_log10_p, color = geneset_color))
g <- g + theme_classic(base_size = 20, base_family = "Helvetica")
g <- g + ggrastr::geom_point_rast(alpha = 0.3, size = 0.8)
g <- g + geom_hline(yintercept = -log10(genomewide_threshold), linetype = "dashed", color = "red")

if (nrow(label_candidates) > 0) {
  g <- g + geom_point(data = label_candidates, size = 2.5, alpha = 0.9, show.legend = FALSE)
  g <- g + geom_text_repel(data = label_candidates, aes(label = LABEL),
    size = 5, max.overlaps = 100, show.legend = FALSE)
}

g <- g + scale_x_continuous(labels = axis_df$CHROM, breaks = axis_df$center)
g <- g + scale_color_manual(values = palette_values)
g <- g + guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
g <- g + theme(legend.title = element_blank())
g <- g + xlab("Chromosome")
g <- g + ylab("-log10(P)")

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = 16, height = 7, dpi = 300)
ggsave(plot = g, filename = paste0(plot_prefix, ".png"), width = 16, height = 7, dpi = 300)
