args <- commandArgs(trailingOnly = TRUE)

gwas_path <- as.character(args[1])
gene_annotation_path <- as.character(args[2])
geneset_dir <- as.character(args[3])
highlight_genesets_raw <- if (length(args) >= 4) as.character(args[4]) else ""
variants_path <- as.character(args[5])
plot_prefix <- as.character(args[6])
flank_bp <- if (length(args) >= 7) as.numeric(args[7]) else 50000
label_p_threshold <- if (length(args) >= 8) as.numeric(args[8]) else 1e-30
genomewide_threshold <- if (length(args) >= 9) as.numeric(args[9]) else 5e-8
data_genesets_raw <- if (length(args) >= 10 && nzchar(args[10])) as.character(args[10]) else ""
spectra_path <- if (length(args) >= 11 && nzchar(args[11])) as.character(args[11]) else ""
gene_map_path <- if (length(args) >= 12 && nzchar(args[12])) as.character(args[12]) else ""
program_k <- if (length(args) >= 13) as.numeric(args[13]) else 60
top_n_program_genes <- if (length(args) >= 14) as.numeric(args[14]) else 100
sampling_trigger_rows <- if (length(args) >= 15) as.numeric(args[15]) else 100000
sampling_base_points <- if (length(args) >= 16) as.numeric(args[16]) else 50000
sampling_fraction <- if (length(args) >= 17) as.numeric(args[17]) else 0.01
sampling_max_points <- if (length(args) >= 18) as.numeric(args[18]) else 300000
sampling_seed <- if (length(args) >= 19) as.numeric(args[19]) else 1

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

parse_name_list <- function(raw) {
  if (is.null(raw) || !nzchar(raw)) {
    return(character())
  }
  values <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  values[nzchar(values)]
}

dedupe_delimited_values <- function(values, sep = ";") {
  vapply(
    values,
    function(value) {
      if (is.na(value) || !nzchar(value)) {
        return("")
      }
      parts <- trimws(strsplit(value, sep, fixed = TRUE)[[1]])
      parts <- parts[nzchar(parts)]
      if (length(parts) == 0) {
        return("")
      }
      paste(unique(parts), collapse = sep)
    },
    character(1)
  )
}

sample_gwas_points <- function(
  gwas_df,
  trigger_rows,
  base_points,
  fraction,
  max_points,
  genomewide_threshold,
  label_p_threshold,
  seed
) {
  total_rows <- nrow(gwas_df)
  target_rows <- total_rows
  sampling_applied <- FALSE

  if (total_rows <= trigger_rows) {
    return(list(
      data = gwas_df,
      total_rows = total_rows,
      target_rows = total_rows,
      sampled_rows = total_rows,
      sampling_applied = sampling_applied
    ))
  }

  target_rows <- base_points + ceiling(total_rows * fraction)
  target_rows <- min(target_rows, max_points)
  target_rows <- min(target_rows, total_rows)
  if (target_rows >= total_rows) {
    return(list(
      data = gwas_df,
      total_rows = total_rows,
      target_rows = total_rows,
      sampled_rows = total_rows,
      sampling_applied = sampling_applied
    ))
  }

  keep_mask <- gwas_df$P <= genomewide_threshold | gwas_df$P <= label_p_threshold
  keep_indices <- which(keep_mask)
  if (length(keep_indices) >= target_rows) {
    message(sprintf(
      "Sampling target %d rows is smaller than mandatory kept rows %d; truncating kept rows by significance.",
      target_rows,
      length(keep_indices)
    ))
    keep_df <- gwas_df[keep_indices, , drop = FALSE]
    keep_df <- keep_df[order(keep_df$P, keep_df$geneset_color == "background", keep_df$POS), , drop = FALSE]
    sampled_df <- head(keep_df, target_rows)
    sampled_df <- sampled_df[order(as.numeric(sampled_df$CHROM), sampled_df$POS), , drop = FALSE]
    return(list(
      data = sampled_df,
      total_rows = total_rows,
      target_rows = target_rows,
      sampled_rows = nrow(sampled_df),
      sampling_applied = TRUE
    ))
  }

  background_indices <- setdiff(seq_len(total_rows), keep_indices)
  set.seed(as.integer(seed))
  sampled_background <- sample(background_indices, target_rows - length(keep_indices), replace = FALSE)
  sampled_indices <- sort(c(keep_indices, sampled_background))
  sampled_df <- gwas_df[sampled_indices, , drop = FALSE]
  sampled_df <- sampled_df[order(as.numeric(sampled_df$CHROM), sampled_df$POS), , drop = FALSE]
  list(
    data = sampled_df,
    total_rows = total_rows,
    target_rows = target_rows,
    sampled_rows = nrow(sampled_df),
    sampling_applied = TRUE
  )
}

build_region_table <- function(annotation_df, geneset_names, geneset_dir, flank_bp) {
  region_rows <- list()
  for (geneset_name in geneset_names) {
    genes <- read_geneset_members(geneset_dir, geneset_name)
    gene_rows <- annotation_df[annotation_df$gene %in% genes, c("gene", "chr_norm", "start", "end"), drop = FALSE]
    if (nrow(gene_rows) == 0) {
      next
    }

    gene_rows <- gene_rows[!duplicated(gene_rows$gene), , drop = FALSE]
    gene_rows$region_start <- pmax(1, gene_rows$start - flank_bp)
    gene_rows$region_end <- gene_rows$end + flank_bp
    gene_rows$geneset_name <- geneset_name
    gene_rows$display_name <- friendly_geneset_name(geneset_name)
    region_rows[[length(region_rows) + 1]] <- gene_rows
  }

  if (length(region_rows) == 0) {
    return(NULL)
  }
  do.call(rbind, region_rows)
}

annotate_nearest_genes <- function(chrom_vec, pos_vec, genes_by_chr) {
  nearest_gene <- rep("", length(chrom_vec))
  distance_to_gene <- rep(NA_real_, length(chrom_vec))

  for (index in seq_along(chrom_vec)) {
    chr_key <- as.character(chrom_vec[index])
    chr_genes <- genes_by_chr[[chr_key]]
    if (is.null(chr_genes) || nrow(chr_genes) == 0) {
      next
    }

    pos_value <- as.numeric(pos_vec[index])
    distances <- ifelse(
      pos_value < chr_genes$start,
      chr_genes$start - pos_value,
      ifelse(pos_value > chr_genes$end, pos_value - chr_genes$end, 0)
    )
    nearest_index <- which.min(distances)
    nearest_gene[index] <- chr_genes$gene[nearest_index]
    distance_to_gene[index] <- distances[nearest_index]
  }

  data.frame(
    nearest_gene = nearest_gene,
    distance_to_gene = distance_to_gene,
    stringsAsFactors = FALSE
  )
}

highlight_genesets <- parse_name_list(highlight_genesets_raw)
data_genesets <- parse_name_list(data_genesets_raw)

gwas <- data.table::fread(gwas_path, data.table = FALSE)
gwas_df <- parse_variant_table(gwas)
coord_data <- add_cumulative_coordinates(gwas_df)
gwas_df <- coord_data$data
axis_df <- coord_data$axis

gene_annotation <- read_gene_annotation_table(gene_annotation_path)
gene_annotation$chr_norm <- sub("^chr", "", as.character(gene_annotation$chr))
gene_annotation <- gene_annotation[
  gene_annotation$chr_norm %in% as.character(seq_len(22)),
  c("gene", "chr_norm", "start", "end"),
  drop = FALSE
]
gene_annotation <- gene_annotation[order(as.numeric(gene_annotation$chr_norm), gene_annotation$start, gene_annotation$end), , drop = FALSE]
genes_by_chr <- split(gene_annotation, gene_annotation$chr_norm)

program_lookup <- character()
if (nzchar(spectra_path) && nzchar(gene_map_path) && file.exists(spectra_path) && file.exists(gene_map_path)) {
  program_lookup <- read_program_membership_lookup(
    spectra_path = spectra_path,
    gene_map_path = gene_map_path,
    k = as.integer(program_k),
    top_n = as.integer(top_n_program_genes)
  )
} else {
  message("Program membership annotation disabled: spectra_path or gene_map_path is unavailable.")
}

region_df <- build_region_table(gene_annotation, highlight_genesets, geneset_dir, flank_bp)
data_region_df <- build_region_table(gene_annotation, data_genesets, geneset_dir, flank_bp)
has_genesets <- !is.null(region_df) && nrow(region_df) > 0
has_data_genesets <- !is.null(data_region_df) && nrow(data_region_df) > 0

gwas_df$geneset <- ""
gwas_df$geneset_color <- "background"

if (has_genesets) {
  for (chr_value in unique(gwas_df$CHROM)) {
    chr_mask <- gwas_df$CHROM == chr_value
    chr_regions <- region_df[region_df$chr_norm == chr_value, , drop = FALSE]
    if (nrow(chr_regions) == 0) {
      next
    }

    chr_regions <- chr_regions[nrow(chr_regions):1, , drop = FALSE]
    chr_pos <- gwas_df$POS[chr_mask]
    chr_geneset_color <- rep("background", sum(chr_mask))
    for (index in seq_len(nrow(chr_regions))) {
      in_region <- chr_pos >= chr_regions$region_start[index] & chr_pos <= chr_regions$region_end[index]
      chr_geneset_color[in_region] <- chr_regions$display_name[index]
    }
    gwas_df$geneset_color[chr_mask] <- chr_geneset_color
  }
}

if (has_data_genesets) {
  for (chr_value in unique(gwas_df$CHROM)) {
    chr_mask <- gwas_df$CHROM == chr_value
    chr_regions <- data_region_df[data_region_df$chr_norm == chr_value, , drop = FALSE]
    if (nrow(chr_regions) == 0) {
      next
    }

    chr_regions <- chr_regions[nrow(chr_regions):1, , drop = FALSE]
    chr_pos <- gwas_df$POS[chr_mask]
    chr_geneset <- rep("", sum(chr_mask))
    for (index in seq_len(nrow(chr_regions))) {
      in_region <- chr_pos >= chr_regions$region_start[index] & chr_pos <= chr_regions$region_end[index]
      chr_geneset[in_region] <- ifelse(
        nchar(chr_geneset[in_region]) == 0,
        chr_regions$display_name[index],
        paste(chr_geneset[in_region], chr_regions$display_name[index], sep = ";")
      )
    }
    gwas_df$geneset[chr_mask] <- chr_geneset
  }
}

gwas_df$geneset <- dedupe_delimited_values(gwas_df$geneset)
gwas_df$geneset[gwas_df$geneset == ""] <- "background"
gwas_df$neg_log10_p <- -log10(pmax(gwas_df$P, .Machine$double.xmin))

sampling_result <- sample_gwas_points(
  gwas_df = gwas_df,
  trigger_rows = as.integer(sampling_trigger_rows),
  base_points = as.integer(sampling_base_points),
  fraction = as.numeric(sampling_fraction),
  max_points = as.integer(sampling_max_points),
  genomewide_threshold = as.numeric(genomewide_threshold),
  label_p_threshold = as.numeric(label_p_threshold),
  seed = as.integer(sampling_seed)
)
gwas_plot_df <- sampling_result$data
if (sampling_result$sampling_applied) {
  message(sprintf(
    "GWAS Manhattan sampling applied: total_rows=%d, target_rows=%d, sampled_rows=%d",
    sampling_result$total_rows,
    sampling_result$target_rows,
    sampling_result$sampled_rows
  ))
}

label_candidates <- gwas_plot_df[
  gwas_plot_df$P <= label_p_threshold & gwas_plot_df$geneset_color != "background",
  c("CHROM", "POS", "SNP", "P", "COORD", "geneset_color", "neg_log10_p"),
  drop = FALSE
]
if (nrow(label_candidates) > 0) {
  label_candidates <- label_candidates[!duplicated(label_candidates$COORD), , drop = FALSE]
  label_info <- annotate_nearest_genes(label_candidates$CHROM, label_candidates$POS, genes_by_chr)
  label_candidates$LABEL <- label_info$nearest_gene
  label_candidates <- label_candidates[label_candidates$LABEL != "", , drop = FALSE]
} else {
  label_candidates <- data.frame(
    CHROM = character(),
    POS = numeric(),
    SNP = character(),
    P = numeric(),
    COORD = numeric(),
    geneset_color = character(),
    neg_log10_p = numeric(),
    LABEL = character(),
    stringsAsFactors = FALSE
  )
}

variants_export <- data.frame(
  chr = gwas_plot_df$CHROM,
  bp = gwas_plot_df$POS,
  snp = gwas_plot_df$SNP,
  p = gwas_plot_df$P,
  logp = gwas_plot_df$neg_log10_p,
  stringsAsFactors = FALSE
)
variants_export$total_input_rows <- sampling_result$total_rows
variants_export$target_sample_rows <- sampling_result$target_rows
variants_export$sampled_output_rows <- sampling_result$sampled_rows
variants_export$sampling_applied <- sampling_result$sampling_applied

hit_mask <- gwas_df$P <= genomewide_threshold
if (any(hit_mask)) {
  hits_df <- gwas_df[hit_mask, c("CHROM", "POS", "SNP", "P", "neg_log10_p", "geneset"), drop = FALSE]
  hit_info <- annotate_nearest_genes(hits_df$CHROM, hits_df$POS, genes_by_chr)
  hit_programs <- unname(program_lookup[hit_info$nearest_gene])
  hit_programs[is.na(hit_programs)] <- ""
  hit_genesets <- hits_df$geneset
  hit_genesets[hit_genesets == "background"] <- ""

  hits_export <- data.frame(
    chr = hits_df$CHROM,
    bp = hits_df$POS,
    snp = hits_df$SNP,
    p = hits_df$P,
    logp = hits_df$neg_log10_p,
    nearest_gene = hit_info$nearest_gene,
    distance_to_gene = hit_info$distance_to_gene,
    program = hit_programs,
    geneset = hit_genesets,
    stringsAsFactors = FALSE
  )
} else {
  hits_export <- data.frame(
    chr = character(),
    bp = numeric(),
    snp = character(),
    p = numeric(),
    logp = numeric(),
    nearest_gene = character(),
    distance_to_gene = numeric(),
    program = character(),
    geneset = character(),
    stringsAsFactors = FALSE
  )
}

hits_path <- gsub("_variants\\.tsv$", "_hits.tsv", variants_path)
if (hits_path == variants_path) {
  hits_path <- gsub("\\.tsv$", "_hits.tsv", variants_path)
}
if (hits_path == variants_path) {
  hits_path <- paste0(variants_path, "_hits.tsv")
}

ensure_parent_dir(variants_path)
ensure_parent_dir(paste0(plot_prefix, ".pdf"))

utils::write.table(variants_export, variants_path, row.names = FALSE, sep = "\t", quote = FALSE)
utils::write.table(hits_export, hits_path, row.names = FALSE, sep = "\t", quote = FALSE)

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
  for (index in seq_along(unique_sets)) {
    palette_values[unique_sets[index]] <- highlight_colors[((index - 1) %% length(highlight_colors)) + 1]
  }
}

g <- ggplot(gwas_plot_df, aes(x = COORD, y = neg_log10_p, color = geneset_color))
g <- g + theme_classic(base_size = 20, base_family = "Helvetica")
if (nrow(gwas_plot_df) > 8000000) {
  g <- g + geom_point(alpha = 0.1, size = 0.3, shape = ".")
} else {
  g <- g + ggrastr::geom_point_rast(alpha = 0.3, size = 0.8, raster.dpi = 200)
}
g <- g + geom_hline(yintercept = -log10(genomewide_threshold), linetype = "dashed", color = "red")

if (nrow(label_candidates) > 0) {
  g <- g + geom_point(data = label_candidates, size = 2.5, alpha = 0.9, show.legend = FALSE)
  g <- g + geom_text_repel(
    data = label_candidates,
    aes(label = LABEL),
    size = 5,
    max.overlaps = 100,
    show.legend = FALSE
  )
}

g <- g + scale_x_continuous(labels = axis_df$CHROM, breaks = axis_df$center)
g <- g + scale_color_manual(values = palette_values)
g <- g + guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
g <- g + theme(legend.title = element_blank())
g <- g + xlab("Chromosome")
g <- g + ylab("-log10(P)")

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = 16, height = 7, dpi = 300)
ggsave(plot = g, filename = paste0(plot_prefix, ".png"), width = 16, height = 7, dpi = 300)
