args <- commandArgs(trailingOnly = TRUE)

posterior_path <- as.character(args[1])
limma_path <- as.character(args[2])
shet_path <- as.character(args[3])
gene_map_path <- as.character(args[4])
trait_label <- as.character(args[5])
table_path <- as.character(args[6])
plot_prefix <- as.character(args[7])
y_limit <- if (length(args) >= 8) as.numeric(args[8]) else 10
render_plot <- if (length(args) >= 9) as.character(args[9]) else "1"

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

options(scipen = 10)
posterior_df <- data.table::fread(posterior_path, data.table = FALSE)
if (!all(c("ensg", "post_mean") %in% names(posterior_df))) {
  stop(sprintf("Posterior file %s must have columns: ensg, post_mean", posterior_path))
}
posterior_df$ensg <- normalize_ensg_ids(posterior_df$ensg)
posterior_df$post_mean[is.infinite(posterior_df$post_mean)] <- NA

limma_raw <- data.table::fread(limma_path, data.table = FALSE)
limma_mat <- as.matrix(limma_raw[, -1, drop = FALSE])
row.names(limma_mat) <- limma_raw[[1]]
target_ids <- row.names(limma_mat)
perturb_ids <- colnames(limma_mat)

shet_df <- utils::read.table(shet_path, header = TRUE, stringsAsFactors = FALSE)
shet_df$ensg <- normalize_ensg_ids(shet_df$ensg)

lookups <- build_gene_id_lookups(gene_map_path)
target_ensg <- lookups$resolve_to_ensg(target_ids)
perturb_ensg <- lookups$resolve_to_ensg(perturb_ids)

shet_sub <- shet_df[shet_df$ensg %in% posterior_df$ensg, ]
n_targets <- nrow(limma_mat)

corr_list <- lapply(seq_len(n_targets), function(i) {
  pb <- limma_mat[i, ]
  pb_df <- data.frame(
    gene = perturb_ids,
    perturb_beta = as.numeric(pb),
    stringsAsFactors = FALSE
  )
  pb_df$ensg <- perturb_ensg

  df <- merge(pb_df, posterior_df[, c("ensg", "post_mean")], by = "ensg", all = FALSE)
  df <- df[!is.na(df$ensg) & !is.na(df$post_mean) & !is.na(df$perturb_beta), ]
  if (!is.na(target_ensg[i])) {
    df <- df[df$ensg != target_ensg[i], ]
  }
  df <- merge(df, shet_sub, by = "ensg", all = FALSE)
  if (nrow(df) < 10) {
    return(NULL)
  }
  if (stats::sd(df$post_mean, na.rm = TRUE) == 0 || stats::sd(df$perturb_beta, na.rm = TRUE) == 0) {
    return(NULL)
  }

  df$post_mean_sc <- scale(df$post_mean)
  df$perturb_beta_sc <- scale(df$perturb_beta)

  fit <- tryCatch(lm(post_mean_sc ~ perturb_beta_sc + shet, data = df), error = function(e) NULL)
  if (is.null(fit)) {
    return(NULL)
  }
  sm <- summary(fit)$coefficients
  if (!"perturb_beta_sc" %in% row.names(sm)) {
    return(NULL)
  }

  data.frame(
    ensg = target_ensg[i],
    P_withShet = sm["perturb_beta_sc", 4],
    beta_withShet = sm["perturb_beta_sc", 1],
    stringsAsFactors = FALSE
  )
})

correlation_df <- do.call(rbind, corr_list[!vapply(corr_list, is.null, logical(1))])
if (is.null(correlation_df) || nrow(correlation_df) == 0) {
  unresolved_targets <- unique(target_ids[is.na(target_ensg)])
  unresolved_perturb <- unique(perturb_ids[is.na(perturb_ensg)])
  stop(sprintf(
    paste(
      "No gene-level QQ points produced valid correlation results",
      "(limma rows=%d, resolved target ENSG=%d, perturb columns=%d, resolved perturb ENSG=%d, posterior genes=%d, shet genes=%d).",
      "Example unresolved target IDs: %s.",
      "Example unresolved perturb IDs: %s."
    ),
    length(target_ids),
    sum(!is.na(target_ensg)),
    length(perturb_ids),
    sum(!is.na(perturb_ensg)),
    length(unique(posterior_df$ensg)),
    length(unique(shet_sub$ensg)),
    if (length(unresolved_targets) > 0) paste(head(unresolved_targets, 5), collapse = ", ") else "none",
    if (length(unresolved_perturb) > 0) paste(head(unresolved_perturb, 5), collapse = ", ") else "none"
  ))
}

correlation_df <- correlation_df[!is.na(correlation_df$ensg) & nzchar(correlation_df$ensg), , drop = FALSE]
if (nrow(correlation_df) == 0) {
  stop("No QQ rows remained after ENSG resolution")
}

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

correlation_df$expected <- signed_expected(n)
correlation_df$observed <- correlation_df$signed_log10_p
correlation_df$observed[correlation_df$observed > y_limit] <- y_limit
correlation_df$observed[correlation_df$observed < -y_limit] <- -y_limit
correlation_df$trait_id <- trait_label
correlation_df$qq_rank <- seq_len(nrow(correlation_df))
correlation_df$tail_side <- ifelse(correlation_df$observed >= 0, "positive", "negative")

ensure_parent_dir(table_path)
utils::write.table(correlation_df, table_path, row.names = FALSE, sep = "\t", quote = FALSE)

if (render_plot %in% c("1", "true", "TRUE", "yes", "YES")) {
  ensure_parent_dir(paste0(plot_prefix, ".pdf"))

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
}
