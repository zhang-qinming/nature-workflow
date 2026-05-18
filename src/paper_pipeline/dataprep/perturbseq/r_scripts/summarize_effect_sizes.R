args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: summarize_effect_sizes.R <limma_dir> <summary_dir>")
}

limma_dir <- args[1]
summary_dir <- args[2]

dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

file_list <- dir(limma_dir)
file_list <- file_list[grep("_KO.txt", file_list)]
if (length(file_list) == 0) {
  stop(paste0("No limma outputs found in ", limma_dir))
}

for (FILE in file_list) {
  data <- read.table(file.path(limma_dir, FILE), header = TRUE, stringsAsFactors = FALSE)
  GENE <- gsub("_KO.txt", "", FILE)
  if (FILE == file_list[1]) {
    beta_mat <- data[, c("GENE", "logFC")]
    colnames(beta_mat)[2] <- GENE
    p_mat <- data[, c("GENE", "P.Value")]
    colnames(p_mat)[2] <- GENE
    adj_p_mat <- data[, c("GENE", "adj.P.Val")]
    colnames(adj_p_mat)[2] <- GENE
  } else {
    beta_mat <- merge(beta_mat, data[, c("GENE", "logFC")], by = "GENE")
    colnames(beta_mat)[ncol(beta_mat)] <- GENE
    p_mat <- merge(p_mat, data[, c("GENE", "P.Value")], by = "GENE")
    colnames(p_mat)[ncol(p_mat)] <- GENE
    adj_p_mat <- merge(adj_p_mat, data[, c("GENE", "adj.P.Val")], by = "GENE")
    colnames(adj_p_mat)[ncol(adj_p_mat)] <- GENE
  }
}

write.table(
  beta_mat,
  file.path(summary_dir, "limma_logFC_sum.txt"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)
write.table(
  p_mat,
  file.path(summary_dir, "limma_P_sum.txt"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)
write.table(
  adj_p_mat,
  file.path(summary_dir, "limma_adjP_sum.txt"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)
