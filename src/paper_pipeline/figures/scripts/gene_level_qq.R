args <- commandArgs(trailingOnly = TRUE)

posterior_path <- as.character(args[1])
limma_path      <- as.character(args[2])
shet_path       <- as.character(args[3])
gene_map_path   <- as.character(args[4])
trait_label     <- as.character(args[5])
table_path      <- as.character(args[6])
plot_prefix     <- as.character(args[7])
y_limit         <- if (length(args) >= 8) as.numeric(args[8]) else 10

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

# ── 读输入 ──
options(scipen = 10)
posterior_df <- data.table::fread(posterior_path, data.table = FALSE)
if (!all(c("ensg", "post_mean") %in% names(posterior_df))) {
  stop(sprintf("Posterior file %s must have columns: ensg, post_mean", posterior_path))
}
posterior_df$post_mean[is.infinite(posterior_df$post_mean)] <- NA

limma_raw <- data.table::fread(limma_path, data.table = FALSE)
limma_mat <- as.matrix(limma_raw[, -1, drop = FALSE])
row.names(limma_mat) <- limma_raw[[1]]
perturb_genes <- colnames(limma_mat)

shet_df <- utils::read.table(shet_path, header = TRUE, stringsAsFactors = FALSE)

gene_map <- read_gene_map(gene_map_path)
ensg_lookup <- stats::setNames(gene_map$ensg, gene_map$gene)

# ── 对每个扰动基因算相关性 ──
shet_sub <- shet_df[shet_df$ensg %in% posterior_df$ensg, ]
n_genes  <- length(perturb_genes)

corr_list <- lapply(seq_len(n_genes), function(i) {
  target_gene <- perturb_genes[i]
  pb <- limma_mat[, i]
  pb_df <- data.frame(
    gene         = names(pb),
    perturb_beta = as.numeric(pb),
    stringsAsFactors = FALSE
  )
  pb_df$ensg <- ensg_lookup[pb_df$gene]

  df <- merge(pb_df, posterior_df[, c("ensg", "post_mean")], by = "ensg", all = FALSE)
  df <- df[!is.na(df$post_mean) & !is.na(df$perturb_beta), ]
  df <- df[df$ensg != ensg_lookup[target_gene], ]
  df <- merge(df, shet_sub, by = "ensg", all = FALSE)
  if (nrow(df) < 10) return(NULL)

  df$post_mean_sc     <- scale(df$post_mean)
  df$perturb_beta_sc  <- scale(df$perturb_beta)

  fit <- lm(post_mean_sc ~ perturb_beta_sc + shet, data = df)
  sm  <- summary(fit)$coefficients

  data.frame(
    ensg          = ensg_lookup[target_gene],
    P_withShet    = sm[2, 4],
    beta_withShet = sm[2, 1],
    stringsAsFactors = FALSE
  )
})

correlation_df <- do.call(rbind, corr_list[!vapply(corr_list, is.null, logical(1))])
if (is.null(correlation_df) || nrow(correlation_df) == 0) {
  stop("No perturb genes produced valid correlation results")
}

# ── 计算 signed log10(P) 和理论分位数 ──
correlation_df$signed_log10_p <- sign(correlation_df$beta_withShet) *
  (-log10(pmax(correlation_df$P_withShet, .Machine$double.xmin)))
correlation_df <- correlation_df[order(correlation_df$signed_log10_p), , drop = FALSE]

n <- nrow(correlation_df)
signed_expected <- function(n) {
  if (n %% 2 == 0) {
    sort(c(-log10(seq(1 / n, 1, by = 2 / n)), log10(seq(1 / n, 1, by = 2 / n))))
  } else {
    sort(c(-log10(seq(1 / n, 1, by = 2 / n)), log10(seq(2 / n, 1, by = 2 / n))))
  }
}

correlation_df$expected  <- signed_expected(n)
correlation_df$observed  <- correlation_df$signed_log10_p
correlation_df$observed[correlation_df$observed >  y_limit] <-  y_limit
correlation_df$observed[correlation_df$observed < -y_limit] <- -y_limit
correlation_df$trait_id  <- trait_label

ensure_parent_dir(table_path)
ensure_parent_dir(paste0(plot_prefix, ".pdf"))
utils::write.table(correlation_df, table_path, row.names = FALSE, sep = "\t", quote = FALSE)

# ── 画图（单 trait QQ） ──
library(ggplot2)

g <- ggplot(correlation_df, aes(x = expected, y = observed))
g <- g + theme_classic(base_size = 18, base_family = "Helvetica")
g <- g + geom_point(alpha = 0.6, size = 2.2, color = "#4C78A8")
g <- g + geom_abline(intercept = 0, slope = 1, linetype = "dashed")
g <- g + xlab("Expected signed -log10(P)")
g <- g + ylab("Observed signed -log10(P)")
g <- g + ggtitle(trait_label)

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = 10, height = 9)
ggsave(plot = g, filename = paste0(plot_prefix, ".png"), width = 10, height = 9, dpi = 300)
