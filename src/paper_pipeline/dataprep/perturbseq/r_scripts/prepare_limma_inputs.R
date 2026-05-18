args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: prepare_limma_inputs.R <metadata_csv> <count_csv> <output_dir> [chunk_size]")
}

options("scipen" = 10)

metadata_path <- args[1]
count_path <- args[2]
output_dir <- args[3]
chunk_size <- if (length(args) >= 4) as.numeric(args[4]) else 50

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

metadata <- read.csv(metadata_path, row.names = 1)

library(data.table)
library(dplyr)

count <- fread(count_path, sep = ",", header = TRUE)
if (!"cell_barcode" %in% colnames(count)) {
  colnames(count)[1] <- "cell_barcode"
}

required_columns <- c("gem_group", "n_genes", "mitopercent", "gene")
missing_columns <- setdiff(required_columns, colnames(metadata))
if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Perturbseq metadata is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

if (!"leiden" %in% colnames(metadata)) {
  metadata$leiden <- NA_character_
}

control <- row.names(metadata %>% filter(gene %in% "non-targeting"))
B <- count[is.element(count$cell_barcode, control), ]
dfB <- cbind(B, metadata[B$cell_barcode, c("gem_group", "n_genes", "mitopercent", "leiden", "gene")])
dfB <- data.frame(dfB, GROUP = "control")
write.table(
  dfB,
  file.path(output_dir, "control.input"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

genes <- unique(sort(metadata$gene))
genes <- genes[!is.element(genes, "non-targeting")]
total_chunks <- ceiling(length(genes) / chunk_size)

for (AA in seq_len(total_chunks)) {
  START <- chunk_size * (AA - 1) + 1
  END <- min(chunk_size * AA, length(genes))
  current_genes <- genes[c(START:END)]

  target <- row.names(metadata[is.element(metadata$gene, current_genes), ])
  A <- count[is.element(count$cell_barcode, target), ]

  dfA <- cbind(A, metadata[A$cell_barcode, c("gem_group", "n_genes", "mitopercent", "leiden", "gene")])
  dfA <- data.frame(dfA, GROUP = "perturb")

  write.table(
    dfA,
    file.path(output_dir, paste0(AA, ".input")),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}
