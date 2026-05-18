args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: limmatrend.R <chunk_index> <input_dir> <output_dir>")
}

options("scipen" = 10)
AA <- as.numeric(args[1])
input_dir <- args[2]
output_dir <- args[3]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

library(dplyr)
library(kSamples)
library(data.table)
library(limma)
library(edgeR)

dfA <- fread(file.path(input_dir, paste0(AA, ".input")), header = TRUE, data.table = FALSE)
dfB <- fread(file.path(input_dir, "control.input"), header = TRUE, data.table = FALSE)
df_all <- rbind(dfA, dfB)

N <- ncol(df_all) - 6

GENES <- unique(df_all$gene)
GENES <- GENES[!is.element(GENES, "non-targeting")]

for (GENE in GENES) {
  out_path <- file.path(output_dir, paste0(GENE, "_KO.txt"))
  if (!file.exists(out_path)) {
    df <- df_all[is.element(df_all$gene, c(GENE, "non-targeting")), ]
    df$gem_group <- factor(as.character(df$gem_group))

    exp <- t(df[, 2:N])
    colnames(exp) <- df[, 1]
    COND <- cbind(model.matrix(~GROUP + gem_group, df), df[, c("n_genes", "mitopercent")])

    dge <- DGEList(exp, group = COND$GROUPperturb)
    dge <- calcNormFactors(dge)
    y <- new("EList")
    y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

    fit <- lmFit(y, design = COND)
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    result <- as.data.frame(topTable(fit, coef = "GROUPperturb", sort.by = "P", n = Inf))

    write.table(
      data.frame(GENE = row.names(result), result),
      out_path,
      row.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
  }
}
