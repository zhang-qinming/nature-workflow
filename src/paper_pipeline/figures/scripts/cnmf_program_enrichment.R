args <- commandArgs(trailingOnly = TRUE)

spectra_path <- as.character(args[1])
gene_map_path <- as.character(args[2])
k <- as.numeric(args[3])
programs_raw <- as.character(args[4])
top_n <- as.numeric(args[5])
geneset_dir <- as.character(args[6])
genesets_raw <- as.character(args[7])
table_prefix <- as.character(args[8])
plot_prefix <- as.character(args[9])

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

genesets <- trimws(strsplit(genesets_raw, ",", fixed = TRUE)[[1]])
genesets <- genesets[nzchar(genesets)]
if (length(genesets) == 0) {
  stop("cnmf_program_enrichment.R requires at least one geneset name to test against")
}

# Read spectra matrix
GEP_raw <- data.table::fread(spectra_path, data.table = FALSE)
if (is.character(GEP_raw[,1])) {
    row.names(GEP_raw) <- GEP_raw[,1]
    GEP_raw <- GEP_raw[,-1]
}
GEP <- t(GEP_raw)
colnames(GEP) <- paste0("P", seq_len(k))

gene_map <- read_gene_map(gene_map_path)
corresp <- gene_map[!duplicated(gene_map$ensg), ]
row.names(corresp) <- corresp$ensg

gene_symbols <- corresp[row.names(GEP), "gene"]
gene_symbols[is.na(gene_symbols)] <- row.names(GEP)[is.na(gene_symbols)]
GEP_df <- data.frame(GENE = gene_symbols, GEP, stringsAsFactors = FALSE)

# Total background size for hypergeometric test
background_genes <- unique(GEP_df$GENE)
total_bg <- length(background_genes)

# Load genesets
geneset_list <- list()
for (gs in genesets) {
  members <- read_geneset_members(geneset_dir, gs)
  # only keep members in our background
  members <- intersect(members, background_genes)
  if (length(members) > 0) {
    geneset_list[[gs]] <- members
  }
}

enrichment_results <- data.frame()

for (prog in programs) {
  tmp <- GEP_df[, c("GENE", prog)]
  colnames(tmp) <- c("GENE", "weight")
  tmp <- tmp[order(tmp$weight, decreasing = TRUE), ]
  prog_top_genes <- head(tmp$GENE, top_n)
  
  for (gs_name in names(geneset_list)) {
    gs_members <- geneset_list[[gs_name]]
    overlap <- intersect(prog_top_genes, gs_members)
    n_overlap <- length(overlap)
    n_gs <- length(gs_members)
    
    # Hypergeometric test
    # phyper(q, m, n, k, lower.tail = FALSE)
    # q = n_overlap - 1
    # m = size of geneset (n_gs)
    # n = size of background not in geneset (total_bg - n_gs)
    # k = sample size (top_n)
    pval <- stats::phyper(n_overlap - 1, n_gs, total_bg - n_gs, top_n, lower.tail = FALSE)
    
    enrichment_results <- rbind(
      enrichment_results,
      data.frame(
        Program = prog,
        Geneset = friendly_geneset_name(gs_name),
        OverlapCount = n_overlap,
        GenesetSize = n_gs,
        P = pval,
        NegLog10P = -log10(max(pval, .Machine$double.xmin)),
        stringsAsFactors = FALSE
      )
    )
  }
}

# Multiple testing correction across all tests
enrichment_results$FDR <- stats::p.adjust(enrichment_results$P, method = "BH")

ensure_parent_dir(paste0(table_prefix, "_enrichment.tsv"))
ensure_parent_dir(paste0(plot_prefix, "_enrichment.pdf"))
utils::write.table(enrichment_results, paste0(table_prefix, "_enrichment.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

library(ggplot2)

enrichment_results$Program <- factor(enrichment_results$Program, levels = programs)

# We can plot a bubble plot (dot plot) where color is -log10P and size is overlap
g <- ggplot(enrichment_results, aes(x = Program, y = Geneset))
g <- g + theme_minimal(base_size = 14, base_family = "Helvetica")
g <- g + geom_point(aes(size = OverlapCount, color = NegLog10P))
g <- g + scale_color_gradient(low = "grey80", high = "#B40426", name = "-log10(P)")
g <- g + scale_size_continuous(name = "Overlap Genes")
g <- g + theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  panel.grid.major = element_line(color = "grey90"),
  panel.border = element_rect(color = "black", fill = NA)
)
g <- g + labs(x = "Program", y = "Gene Set")

plot_width <- max(8, length(programs) * 0.4 + 4)
plot_height <- max(5, length(genesets) * 0.4 + 2)

ggsave(plot = g, filename = paste0(plot_prefix, "_enrichment.pdf"), width = plot_width, height = plot_height)
ggsave(plot = g, filename = paste0(plot_prefix, "_enrichment.png"), width = plot_width, height = plot_height, dpi = 300)
