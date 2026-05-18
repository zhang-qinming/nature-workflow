args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: cnmf_regulatory_effects.R <metadata.csv> <usages.txt> <program_index> <output.txt>")
}

options("scipen" = 10)

metadata_path <- args[1]
usages_path <- args[2]
program_index <- as.numeric(args[3])
output_path <- args[4]
column_index <- program_index + 1

metadata <- read.csv(metadata_path, row.names = 1)

library(data.table)
count <- fread(usages_path, header = TRUE, data.table = FALSE)

GENES <- unique(sort(metadata$gene))
GENES <- GENES[!is.element(GENES, "non-targeting")]

library(dplyr)
control <- row.names(metadata %>% filter(gene %in% "non-targeting"))

summary <- data.frame()

B <- count[is.element(count$cell_barcode, control), c(1, column_index)]
dfB <- cbind(B, metadata[B$cell_barcode, c("gem_group", "n_genes", "mitopercent")])
dfB <- data.frame(dfB, GROUP = "control")

A <- count[, c(1, column_index)]
tmp <- intersect(unique(A$cell_barcode), row.names(metadata))
row.names(A) <- A$cell_barcode
A <- A[tmp, ]
A <- cbind(A, metadata[tmp, c("gem_group", "n_genes", "mitopercent")])

for (GENE in GENES) {
  target <- row.names(metadata %>% filter(gene %in% GENE))
  if (length(target) > 10) {
    dfA <- A[is.element(A$cell_barcode, target), ]
    dfA <- data.frame(dfA, GROUP = "perturb")

    if (length(which(is.element(dfB$gem_group, unique(dfA$gem_group)))) > 0) {
      df <- rbind(dfA, dfB)
      df$gem_group <- factor(as.character(df$gem_group))

      colnames(df)[2] <- "exp"
      df$exp <- scale(df$exp)

      m.lm <- lm(exp ~ GROUP + gem_group + n_genes + mitopercent, data = df)
      lm_es <- summary(m.lm)$coefficients[2, 1]
      lm_p <- summary(m.lm)$coefficients[2, 4]
      hoge <- data.frame(GENE = GENE, lm_es = lm_es, lm_p = lm_p)
      summary <- rbind(summary, hoge)
    }
  }
}

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
write.table(summary, output_path, row.names = FALSE, sep = "\t", quote = FALSE, append = FALSE)
