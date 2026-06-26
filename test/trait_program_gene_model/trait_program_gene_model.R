args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 14) {
  stop(paste(
    "usage: trait_program_gene_model.R <regulation_dir> <spectra_path> <gene_map_path>",
    "<posterior_path> <trait_id> <k> <table_prefix> <plot_prefix> <shet_path>",
    "<program_n> <regulator_n> <loading_top_n> <hit_abs_gamma_threshold>",
    "<regulator_fdr_threshold> [random_iterations] [permutation_iterations] [max_genes_per_side] [render_plot]",
    "[program_fdr_threshold] [regulator_max_n]"
  ))
}

regulation_dir <- as.character(args[1])
spectra_path <- as.character(args[2])
gene_map_path <- as.character(args[3])
posterior_path <- as.character(args[4])
trait_id <- as.character(args[5])
k <- as.integer(args[6])
table_prefix <- as.character(args[7])
plot_prefix <- as.character(args[8])
shet_path <- as.character(args[9])
program_n_arg <- as.character(args[10])
regulator_n_arg <- as.character(args[11])
loading_top_n <- as.integer(args[12])
hit_abs_gamma_threshold <- as.numeric(args[13])
regulator_fdr_threshold <- as.numeric(args[14])
random_iterations <- if (length(args) >= 15) as.integer(args[15]) else 10000L
permutation_iterations <- if (length(args) >= 16) as.integer(args[16]) else 0L
max_genes_per_side <- if (length(args) >= 17) as.integer(args[17]) else 8L
render_plot <- if (length(args) >= 18) as.character(args[18]) else "1"
program_fdr_threshold <- if (length(args) >= 19) as.numeric(args[19]) else 0.05
regulator_max_n <- if (length(args) >= 20) as.integer(args[20]) else 8L
plot_enabled <- render_plot %in% c("1", "true", "TRUE", "yes", "YES")

options(scipen = 10)

parse_selection_count <- function(value, name) {
  value <- trimws(as.character(value))
  if (tolower(value) %in% c("auto", "automatic", "model", "bic")) {
    return(list(mode = "auto", n = NA_integer_))
  }
  n <- suppressWarnings(as.integer(value))
  if (is.na(n) || n < 0) {
    stop(sprintf("%s must be a non-negative integer or auto: %s", name, value))
  }
  list(mode = "fixed", n = n)
}

program_selection <- parse_selection_count(program_n_arg, "program_n")
regulator_selection <- parse_selection_count(regulator_n_arg, "regulator_n")
if (is.na(program_fdr_threshold) || program_fdr_threshold < 0 || program_fdr_threshold > 1) {
  stop(sprintf("program_fdr_threshold must be in [0, 1]: %s", program_fdr_threshold))
}
if (is.na(regulator_max_n) || regulator_max_n < 0) {
  stop(sprintf("regulator_max_n must be a non-negative integer: %s", regulator_max_n))
}

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

normalize_ensg_ids <- function(ids) {
  ids <- trimws(as.character(ids))
  ids[is.na(ids)] <- ""
  sub("\\..*$", "", ids)
}

normalize_program_ids <- function(programs) {
  raw <- trimws(as.character(programs))
  raw[is.na(raw)] <- ""
  stripped <- sub("^P", "", raw, ignore.case = TRUE)
  numeric_program <- grepl("^[0-9]+$", stripped)
  raw[numeric_program] <- paste0("P", as.integer(stripped[numeric_program]))
  raw
}

sign_to_label <- function(values) {
  values <- as.numeric(values)
  labels <- rep("zero", length(values))
  labels[is.na(values)] <- ""
  labels[!is.na(values) & values > 0] <- "positive"
  labels[!is.na(values) & values < 0] <- "negative"
  labels
}

empty_prediction_table <- function() {
  data.frame(
    trait_id = character(),
    Program = character(),
    side = character(),
    gene = character(),
    ensg = character(),
    post_mean = numeric(),
    abs_gamma = numeric(),
    gamma_sign = character(),
    membership_score = numeric(),
    rank_within_side = numeric(),
    predicted_sign = character(),
    post_mean_sign = character(),
    is_concordant = logical(),
    is_discordant = logical(),
    display_bucket = character(),
    display_bucket_label = character(),
    display_column = character(),
    display_column_rank = numeric(),
    program_label = character(),
    gene_label = character(),
    display_rank = numeric(),
    x = numeric(),
    y = numeric(),
    panel_row = numeric(),
    y_global = numeric(),
    has_overlap = logical(),
    empty_reason = character(),
    predicted_effect = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_regulator_effect_table <- function() {
  data.frame(
    gene = character(),
    ensg = character(),
    Program = character(),
    beta = numeric(),
    p = numeric(),
    fdr = numeric(),
    weighted_effect_program = numeric(),
    weighted_effect = numeric(),
    support_count = integer(),
    membership_score = numeric(),
    stringsAsFactors = FALSE
  )
}

read_gene_map <- function(path) {
  gene_map <- utils::read.table(path, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(gene_map) < 2) {
    stop(sprintf("gene map must have at least two columns: %s", path))
  }
  gene_map <- gene_map[, 1:2]
  colnames(gene_map) <- c("ensg", "gene")
  gene_map$ensg <- normalize_ensg_ids(gene_map$ensg)
  gene_map$gene <- trimws(as.character(gene_map$gene))
  gene_map <- gene_map[nzchar(gene_map$ensg) & nzchar(gene_map$gene), , drop = FALSE]
  gene_map <- gene_map[!duplicated(gene_map$ensg), , drop = FALSE]
  gene_map <- gene_map[!duplicated(gene_map$gene), , drop = FALSE]
  gene_map
}

read_spectra_matrix <- function(path, gene_map_path, k) {
  raw <- data.table::fread(path, data.table = FALSE)
  if (ncol(raw) < 1) {
    stop(sprintf("cNMF spectra is empty: %s", path))
  }
  if (is.character(raw[[1]])) {
    row.names(raw) <- raw[[1]]
    raw <- raw[, -1, drop = FALSE]
  }
  mat <- t(as.matrix(raw))
  storage.mode(mat) <- "double"
  actual_k <- min(as.integer(k), ncol(mat))
  mat <- mat[, seq_len(actual_k), drop = FALSE]
  colnames(mat) <- paste0("P", seq_len(actual_k))

  gene_map <- read_gene_map(gene_map_path)
  gene_lookup <- stats::setNames(gene_map$gene, gene_map$ensg)
  ensg <- normalize_ensg_ids(row.names(mat))
  gene <- unname(gene_lookup[ensg])
  gene[is.na(gene) | !nzchar(gene)] <- ensg[is.na(gene) | !nzchar(gene)]

  list(mat = mat, ensg = ensg, gene = gene, gene_map = gene_map)
}

read_regulation_matrices <- function(regulation_dir, k, gene_map) {
  ensg_lookup <- stats::setNames(gene_map$ensg, gene_map$gene)
  beta_df <- data.frame()
  p_df <- data.frame()

  for (i in seq_len(k)) {
    path <- file.path(regulation_dir, paste0("K", k, "_program", i, "_perturb_effects.txt"))
    tmp <- data.table::fread(path, data.table = FALSE)
    if (ncol(tmp) < 3) {
      stop(sprintf("regulation file must have at least three columns: %s", path))
    }
    colnames(tmp)[1:3] <- c("gene", "beta", "p")
    tmp <- tmp[, c("gene", "beta", "p"), drop = FALSE]
    tmp$gene <- trimws(as.character(tmp$gene))
    tmp <- tmp[nzchar(tmp$gene) & !duplicated(tmp$gene), , drop = FALSE]
    tmp_beta <- tmp[, c("gene", "beta"), drop = FALSE]
    tmp_p <- tmp[, c("gene", "p"), drop = FALSE]
    colnames(tmp_beta) <- c("gene", paste0("P", i))
    colnames(tmp_p) <- c("gene", paste0("P", i))
    if (i == 1) {
      beta_df <- tmp_beta
      p_df <- tmp_p
    } else {
      beta_df <- merge(beta_df, tmp_beta, by = "gene", all = FALSE)
      p_df <- merge(p_df, tmp_p, by = "gene", all = FALSE)
    }
  }

  beta_df$ensg <- unname(ensg_lookup[as.character(beta_df$gene)])
  p_df$ensg <- unname(ensg_lookup[as.character(p_df$gene)])
  beta_df <- beta_df[!is.na(beta_df$ensg) & nzchar(beta_df$ensg), , drop = FALSE]
  p_df <- p_df[!is.na(p_df$ensg) & nzchar(p_df$ensg), , drop = FALSE]
  list(beta = beta_df, p = p_df)
}

rank_program_burden <- function(posterior_df, spectra, shet_df, top_n, random_iterations) {
  lof_df <- posterior_df[, c("ensg", "gene", "post_mean"), drop = FALSE]
  rows <- list()
  for (program in colnames(spectra$mat)) {
    tmp <- data.frame(
      ensg = spectra$ensg,
      gene = spectra$gene,
      GEPscore = as.numeric(spectra$mat[, program]),
      stringsAsFactors = FALSE
    )
    tmp <- tmp[order(tmp$GEPscore, decreasing = TRUE, na.last = NA), , drop = FALSE]
    tmp <- tmp[!duplicated(tmp$ensg) & nzchar(tmp$ensg), , drop = FALSE]
    df <- merge(tmp, lof_df, by = c("ensg", "gene"), all = FALSE)
    df <- df[!is.na(df$post_mean), , drop = FALSE]
    if (nrow(df) == 0) {
      next
    }
    top_df <- head(df[order(df$GEPscore, decreasing = TRUE), , drop = FALSE], top_n)
    if (nrow(top_df) == 0) {
      next
    }
    mean_gamma_top <- mean(top_df$post_mean, na.rm = TRUE)
    shet_top <- shet_df[shet_df$ensg %in% top_df$ensg, , drop = FALSE]
    shet_all <- shet_df[shet_df$ensg %in% df$ensg, , drop = FALSE]
    if (nrow(shet_top) == 0 || nrow(shet_all) == 0 || !"shet_BIN" %in% names(shet_all)) {
      next
    }
    bin_counts <- table(shet_top$shet_BIN)
    random_means <- numeric()
    for (seed in seq_len(random_iterations)) {
      set.seed(seed)
      sampled <- character()
      valid <- TRUE
      for (bin_name in names(bin_counts)) {
        candidates <- shet_all$ensg[shet_all$shet_BIN == bin_name]
        if (length(candidates) < bin_counts[[bin_name]]) {
          valid <- FALSE
          break
        }
        sampled <- c(sampled, sample(candidates, bin_counts[[bin_name]]))
      }
      if (valid) {
        random_means <- c(random_means, mean(df$post_mean[df$ensg %in% sampled], na.rm = TRUE))
      }
    }
    if (length(random_means) == 0) {
      next
    }
    rank_value <- rank(c(mean_gamma_top, random_means))[1]
    p1 <- (rank_value / length(random_means)) * 2
    p2 <- ((length(random_means) + 2 - rank_value) / length(random_means)) * 2
    rows[[length(rows) + 1]] <- data.frame(
      Program = program,
      P = min(1, min(p1, p2)),
      meanG = mean_gamma_top - mean(random_means),
      MEANgamma_top = mean_gamma_top,
      shet_adjusted_random_mean = mean(random_means),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) {
    return(data.frame(
      Program = character(),
      P = numeric(),
      q_value = numeric(),
      meanG = numeric(),
      MEANgamma_top = numeric(),
      shet_adjusted_random_mean = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  out$q_value <- stats::p.adjust(out$P, method = "BH")
  out <- out[order(abs(out$meanG), decreasing = TRUE), , drop = FALSE]
  out <- out[order(out$q_value, out$P), , drop = FALSE]
  out
}

select_programs <- function(program_rank, program_selection, fdr_threshold) {
  if (program_selection$mode == "auto") {
    selected <- program_rank$Program[
      !is.na(program_rank$q_value) &
        program_rank$q_value <= fdr_threshold
    ]
    return(selected)
  }
  head(program_rank$Program, program_selection$n)
}

extract_regsubsets_programs <- function(fit_summary, row_index) {
  fit_which <- fit_summary$which
  selected <- colnames(fit_which)[fit_which[row_index, ]]
  selected <- selected[!selected %in% c("(Intercept)", "shet")]
  normalize_program_ids(selected)
}

select_regulator_model <- function(posterior_df, beta_df, shet_df, regulator_selection, regulator_max_n) {
  model_df <- merge(posterior_df[, c("gene", "ensg", "post_mean"), drop = FALSE], beta_df, by = c("gene", "ensg"), all = FALSE)
  model_df <- merge(model_df, shet_df[, c("ensg", "shet"), drop = FALSE], by = "ensg", all = FALSE)
  program_cols <- grep("^P[0-9]+$", names(model_df), value = TRUE)
  model_df <- model_df[, c("post_mean", program_cols, "shet"), drop = FALSE]
  for (col in names(model_df)) {
    finite_values <- model_df[[col]][is.finite(model_df[[col]])]
    if (length(finite_values) > 0) {
      model_df[[col]][is.infinite(model_df[[col]])] <- max(finite_values)
    }
  }

  suppressPackageStartupMessages(library(leaps))
  predictor_n <- length(program_cols) + 1L
  if (regulator_selection$mode == "auto") {
    nvmax <- min(max(1L, regulator_max_n + 1L), predictor_n)
    fit <- leaps::regsubsets(
      post_mean ~ .,
      data = model_df,
      nbest = 1,
      nvmax = nvmax,
      really.big = TRUE
    )
    fit_summary <- summary(fit)
    best_row <- which.min(fit_summary$bic)
    selected <- extract_regsubsets_programs(fit_summary, best_row)
  } else {
    regulator_n <- regulator_selection$n
    if (regulator_n == 0) {
      selected <- character()
    } else {
      fit <- leaps::regsubsets(
        post_mean ~ .,
        data = model_df,
        nbest = 1,
        nvmax = min(regulator_n + 1L, predictor_n),
        really.big = TRUE
      )
      fit_summary <- summary(fit)
      selected <- extract_regsubsets_programs(fit_summary, nrow(fit_summary$which))

      if (length(selected) != regulator_n) {
        fit <- leaps::regsubsets(
          post_mean ~ .,
          data = model_df,
          nbest = 1,
          nvmax = min(regulator_n, predictor_n),
          really.big = TRUE
        )
        fit_summary <- summary(fit)
        selected <- extract_regsubsets_programs(fit_summary, nrow(fit_summary$which))
      }
    }
  }

  lm_df <- model_df[, c("post_mean", selected, "shet"), drop = FALSE]
  for (col in names(lm_df)) {
    lm_df[[col]] <- as.numeric(scale(lm_df[[col]]))
  }
  lm_fit <- stats::lm(post_mean ~ ., data = lm_df)
  coef_df <- as.data.frame(summary(lm_fit)$coefficients)
  coef_df$term <- row.names(coef_df)
  colnames(coef_df)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
  coef_df <- coef_df[coef_df$term %in% selected, c("term", "estimate", "std_error", "t_value", "p_value"), drop = FALSE]
  colnames(coef_df)[1] <- "Program"

  return(list(selected = selected, coefficients = coef_df))
}

top_loading_membership <- function(spectra, selected_programs, top_n) {
  rows <- list()
  for (program in selected_programs) {
    tmp <- data.frame(
      Program = program,
      ensg = spectra$ensg,
      gene = spectra$gene,
      weight = as.numeric(spectra$mat[, program]),
      stringsAsFactors = FALSE
    )
    tmp <- tmp[order(tmp$weight, decreasing = TRUE, na.last = NA), , drop = FALSE]
    tmp <- tmp[!duplicated(tmp$ensg) & nzchar(tmp$ensg), , drop = FALSE]
    tmp <- head(tmp, top_n)
    tmp$rank_within_side <- seq_len(nrow(tmp))
    rows[[length(rows) + 1]] <- tmp
  }
  if (length(rows) == 0) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

build_gene_predictions <- function(
  posterior_df,
  spectra,
  beta_df,
  p_df,
  program_rank,
  program_selected,
  regulator_selected,
  regulator_coef,
  top_n,
  hit_threshold,
  regulator_fdr_threshold,
  max_genes_per_side
) {
  trait_hits <- posterior_df[!is.na(posterior_df$post_mean) & abs(posterior_df$post_mean) > hit_threshold, , drop = FALSE]
  program_sign <- stats::setNames(sign(program_rank$meanG), program_rank$Program)
  coef_lookup <- stats::setNames(regulator_coef$estimate, regulator_coef$Program)

  loading <- top_loading_membership(spectra, program_selected, top_n)
  loading_hits <- merge(
    loading,
    trait_hits[, c("ensg", "gene", "post_mean"), drop = FALSE],
    by = c("ensg", "gene"),
    all = FALSE
  )
  if (nrow(loading_hits) > 0) {
    loading_hits$side <- "program_loading"
    loading_hits$membership_score <- loading_hits$weight
    loading_hits$predicted_effect <- unname(program_sign[loading_hits$Program])
    loading_hits$rank_effect <- ave(-abs(loading_hits$post_mean), loading_hits$Program, FUN = rank, ties.method = "first")
    loading_hits <- loading_hits[order(loading_hits$Program, loading_hits$rank_effect, loading_hits$rank_within_side), , drop = FALSE]
    loading_hits <- do.call(rbind, lapply(split(loading_hits, loading_hits$Program), function(df) head(df, max_genes_per_side)))
  }

  beta_long <- data.frame()
  p_long <- data.frame()
  if (length(regulator_selected) > 0) {
    for (program in regulator_selected) {
      beta_long <- rbind(beta_long, data.frame(gene = beta_df$gene, ensg = beta_df$ensg, Program = program, beta = beta_df[[program]]))
      p_long <- rbind(p_long, data.frame(gene = p_df$gene, ensg = p_df$ensg, Program = program, p = p_df[[program]]))
    }
    reg_long <- merge(beta_long, p_long, by = c("gene", "ensg", "Program"), all = FALSE)
  } else {
    reg_long <- data.frame(gene = character(), ensg = character(), Program = character(), beta = numeric(), p = numeric())
  }
  if (nrow(reg_long) > 0) {
    reg_long$fdr <- ave(reg_long$p, reg_long$Program, FUN = function(x) stats::p.adjust(x, method = "BH"))
    reg_long <- reg_long[reg_long$fdr <= regulator_fdr_threshold, , drop = FALSE]
    reg_long$weighted_effect <- unname(coef_lookup[reg_long$Program]) * reg_long$beta
  }

  regulator_gene_effect <- empty_regulator_effect_table()
  if (nrow(reg_long) > 0) {
    total_effect <- aggregate(
      weighted_effect ~ gene + ensg,
      data = reg_long,
      FUN = sum
    )
    support_count <- aggregate(Program ~ gene + ensg, data = reg_long, FUN = length)
    colnames(support_count)[3] <- "support_count"
    regulator_gene_effect <- merge(reg_long, total_effect, by = c("gene", "ensg"), all = FALSE, suffixes = c("_program", ""))
    regulator_gene_effect <- merge(regulator_gene_effect, support_count, by = c("gene", "ensg"), all = FALSE)
    regulator_gene_effect$membership_score <- abs(regulator_gene_effect$weighted_effect_program)
  }

  regulator_hits <- merge(
    regulator_gene_effect,
    trait_hits[, c("ensg", "gene", "post_mean"), drop = FALSE],
    by = c("ensg", "gene"),
    all = FALSE
  )
  if (nrow(regulator_hits) > 0) {
    regulator_hits$side <- "regulator"
    regulator_hits$predicted_effect <- regulator_hits$weighted_effect
    regulator_hits$rank_within_side <- rank(-abs(regulator_hits$weighted_effect), ties.method = "first")
    regulator_hits <- regulator_hits[order(regulator_hits$Program, -abs(regulator_hits$post_mean), regulator_hits$rank_within_side), , drop = FALSE]
    regulator_hits <- do.call(rbind, lapply(split(regulator_hits, regulator_hits$Program), function(df) head(df, max_genes_per_side)))
  }

  common_cols <- c(
    "Program", "side", "gene", "ensg", "post_mean", "membership_score",
    "rank_within_side", "predicted_effect"
  )
  parts <- list()
  if (exists("loading_hits") && nrow(loading_hits) > 0) {
    parts[[length(parts) + 1]] <- loading_hits[, common_cols]
  }
  if (exists("regulator_hits") && nrow(regulator_hits) > 0) {
    parts[[length(parts) + 1]] <- regulator_hits[, common_cols]
  }
  if (length(parts) == 0) {
    return(empty_prediction_table())
  }

  out <- do.call(rbind, parts)
  out$trait_id <- trait_id
  out$abs_gamma <- abs(out$post_mean)
  out$gamma_sign <- sign_to_label(out$post_mean)
  out$predicted_sign <- sign_to_label(out$predicted_effect)
  out$post_mean_sign <- sign_to_label(out$post_mean)
  out$is_concordant <- sign(out$predicted_effect) != 0 & sign(out$post_mean) != 0 & sign(out$predicted_effect) == sign(out$post_mean)
  out$is_discordant <- sign(out$predicted_effect) != 0 & sign(out$post_mean) != 0 & sign(out$predicted_effect) == -sign(out$post_mean)
  out$display_bucket <- ifelse(out$side == "program_loading", "program_genes", "regulator_genes")
  out$display_bucket_label <- ifelse(out$side == "program_loading", "Program genes", "Regulator genes")
  out$display_column <- ifelse(out$is_concordant, "left", "right")
  out$display_column_rank <- ave(-out$abs_gamma, out$side, FUN = rank, ties.method = "first")
  out$program_label <- out$Program
  out$gene_label <- ifelse(out$is_discordant, paste0("(", out$gene, ")"), out$gene)
  out$display_rank <- ave(-out$abs_gamma, out$Program, out$side, FUN = rank, ties.method = "first")
  out$x <- ifelse(out$side == "program_loading", -1, 1)
  out$y <- out$display_rank
  out$panel_row <- match(out$Program, unique(out$Program))
  out$y_global <- out$y
  out$has_overlap <- TRUE
  out$empty_reason <- ""
  out[, c(
    "trait_id", "Program", "side", "gene", "ensg", "post_mean", "abs_gamma", "gamma_sign",
    "membership_score", "rank_within_side", "predicted_sign", "post_mean_sign",
    "is_concordant", "is_discordant", "display_bucket", "display_bucket_label",
    "display_column", "display_column_rank", "program_label", "gene_label", "display_rank",
    "x", "y", "panel_row", "y_global", "has_overlap", "empty_reason", "predicted_effect"
  )]
}

run_model <- function(posterior_df, do_permutation = FALSE, seed = NA_integer_) {
  df <- posterior_df
  if (do_permutation) {
    set.seed(seed)
    df$post_mean <- sample(df$post_mean, nrow(df), replace = FALSE)
  }
  program_rank <- rank_program_burden(df, spectra, shet_df, loading_top_n, random_iterations)
  if (nrow(program_rank) == 0) {
    stop("No program-burden rankings were produced; check posterior, spectra, gene map, and shet overlap")
  }
  program_selected <- select_programs(program_rank, program_selection, program_fdr_threshold)
  reg_model <- select_regulator_model(df, reg_mats$beta, shet_df, regulator_selection, regulator_max_n)
  predictions <- build_gene_predictions(
    posterior_df = df,
    spectra = spectra,
    beta_df = reg_mats$beta,
    p_df = reg_mats$p,
    program_rank = program_rank,
    program_selected = program_selected,
    regulator_selected = reg_model$selected,
    regulator_coef = reg_model$coefficients,
    top_n = loading_top_n,
    hit_threshold = hit_abs_gamma_threshold,
    regulator_fdr_threshold = regulator_fdr_threshold,
    max_genes_per_side = max_genes_per_side
  )
  list(
    program_rank = program_rank,
    program_selected = program_selected,
    regulator_selected = reg_model$selected,
    regulator_coefficients = reg_model$coefficients,
    predictions = predictions
  )
}

spectra <- read_spectra_matrix(spectra_path, gene_map_path, k)
reg_mats <- read_regulation_matrices(regulation_dir, k, spectra$gene_map)

posterior_df <- data.table::fread(posterior_path, data.table = FALSE)
if (!all(c("ensg", "post_mean") %in% names(posterior_df))) {
  stop(sprintf("posterior file must have ensg and post_mean columns: %s", posterior_path))
}
posterior_df$ensg <- normalize_ensg_ids(posterior_df$ensg)
posterior_df$post_mean[is.infinite(posterior_df$post_mean)] <- NA
gene_lookup <- stats::setNames(spectra$gene_map$gene, spectra$gene_map$ensg)
posterior_df$gene <- unname(gene_lookup[posterior_df$ensg])
posterior_df$gene[is.na(posterior_df$gene) | !nzchar(posterior_df$gene)] <- posterior_df$ensg[
  is.na(posterior_df$gene) | !nzchar(posterior_df$gene)
]
posterior_df <- posterior_df[!is.na(posterior_df$post_mean), c("ensg", "gene", "post_mean"), drop = FALSE]
posterior_df <- posterior_df[nzchar(posterior_df$gene) & !duplicated(posterior_df$gene), , drop = FALSE]

shet_df <- utils::read.table(shet_path, header = TRUE, stringsAsFactors = FALSE)
shet_df$ensg <- normalize_ensg_ids(shet_df$ensg)
if (!"shet_BIN" %in% names(shet_df)) {
  stop(sprintf("shet file must contain shet_BIN: %s", shet_path))
}
if (!"shet" %in% names(shet_df)) {
  shet_df$shet <- as.numeric(shet_df$shet_BIN)
}
shet_df <- shet_df[shet_df$ensg %in% posterior_df$ensg, , drop = FALSE]

model <- run_model(posterior_df)

program_rows <- model$program_rank
program_rows$selected_by_program <- program_rows$Program %in% model$program_selected
program_rows$selected_by_regulator <- program_rows$Program %in% model$regulator_selected
program_rows$trait_id <- trait_id
program_rows$program_trait_sign <- sign_to_label(program_rows$meanG)
program_rows$program_selection_mode <- program_selection$mode
program_rows$program_fdr_threshold <- ifelse(program_selection$mode == "auto", program_fdr_threshold, NA_real_)
program_rows$regulator_selection_mode <- regulator_selection$mode
program_rows$regulator_max_n <- ifelse(regulator_selection$mode == "auto", regulator_max_n, NA_integer_)

coef_rows <- model$regulator_coefficients
coef_rows$trait_id <- rep(trait_id, nrow(coef_rows))
coef_rows$selected_by_regulator <- rep(TRUE, nrow(coef_rows))
coef_rows$regulator_selection_mode <- rep(regulator_selection$mode, nrow(coef_rows))
coef_rows$regulator_max_n <- rep(ifelse(regulator_selection$mode == "auto", regulator_max_n, NA_integer_), nrow(coef_rows))

selected_program_ids <- unique(c(model$program_selected, model$regulator_selected))
if (length(selected_program_ids) > 0) {
  selected_rows <- data.frame(
    trait_id = trait_id,
    Program = selected_program_ids,
    selected_by_program = selected_program_ids %in% model$program_selected,
    selected_by_regulator = selected_program_ids %in% model$regulator_selected,
    stringsAsFactors = FALSE
  )
} else {
  selected_rows <- data.frame(
    trait_id = character(),
    Program = character(),
    selected_by_program = logical(),
    selected_by_regulator = logical(),
    stringsAsFactors = FALSE
  )
}
selected_rows <- merge(
  selected_rows,
  program_rows[, c("Program", "P", "q_value", "meanG", "MEANgamma_top", "shet_adjusted_random_mean", "program_trait_sign")],
  by = "Program",
  all.x = TRUE
)
selected_rows <- merge(
  selected_rows,
  coef_rows[, c("Program", "estimate", "p_value")],
  by = "Program",
  all.x = TRUE
)
colnames(selected_rows)[colnames(selected_rows) == "estimate"] <- "regulator_model_coef"
colnames(selected_rows)[colnames(selected_rows) == "p_value"] <- "regulator_model_p"
selected_rows$program_selection_mode <- rep(program_selection$mode, nrow(selected_rows))
selected_rows$program_fdr_threshold <- rep(ifelse(program_selection$mode == "auto", program_fdr_threshold, NA_real_), nrow(selected_rows))
selected_rows$regulator_selection_mode <- rep(regulator_selection$mode, nrow(selected_rows))
selected_rows$regulator_max_n <- rep(ifelse(regulator_selection$mode == "auto", regulator_max_n, NA_integer_), nrow(selected_rows))
selected_rows$loading_gene_count <- integer(nrow(selected_rows))
selected_rows$regulator_gene_count <- integer(nrow(selected_rows))

predictions <- model$predictions
graph_long <- predictions[predictions$is_concordant, , drop = FALSE]
if (nrow(graph_long) > 0) {
  loading_counts <- table(graph_long$Program[graph_long$side == "program_loading"])
  regulator_counts <- table(graph_long$Program[graph_long$side == "regulator"])
  selected_rows$loading_gene_count <- as.integer(loading_counts[selected_rows$Program])
  selected_rows$regulator_gene_count <- as.integer(regulator_counts[selected_rows$Program])
  selected_rows$loading_gene_count[is.na(selected_rows$loading_gene_count)] <- 0L
  selected_rows$regulator_gene_count[is.na(selected_rows$regulator_gene_count)] <- 0L
}
selected_rows$program_label <- sprintf(
  "%s  L:%d  R:%d",
  selected_rows$Program,
  selected_rows$loading_gene_count,
  selected_rows$regulator_gene_count
)
selected_rows$has_overlap <- rep(nrow(graph_long) > 0, nrow(selected_rows))
selected_rows$empty_reason <- ifelse(selected_rows$has_overlap, "", "no_concordant_trait_program_gene_overlap")

permutation_rows <- data.frame()
if (permutation_iterations > 0) {
  observed_concordant <- sum(predictions$is_concordant, na.rm = TRUE)
  observed_discordant <- sum(predictions$is_discordant, na.rm = TRUE)
  for (seed in seq_len(permutation_iterations)) {
    perm <- run_model(posterior_df, do_permutation = TRUE, seed = seed)
    perm_pred <- perm$predictions
    permutation_rows <- rbind(
      permutation_rows,
      data.frame(
        trait_id = trait_id,
        seed = seed,
        concordant_n = sum(perm_pred$is_concordant, na.rm = TRUE),
        discordant_n = sum(perm_pred$is_discordant, na.rm = TRUE),
        observed_concordant_n = observed_concordant,
        observed_discordant_n = observed_discordant,
        stringsAsFactors = FALSE
      )
    )
  }
}

ensure_parent_dir(paste0(table_prefix, "_long.tsv"))
utils::write.table(graph_long, paste0(table_prefix, "_long.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
utils::write.table(selected_rows, paste0(table_prefix, "_programs.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
utils::write.table(predictions, paste0(table_prefix, "_gene_predictions.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
utils::write.table(program_rows, paste0(table_prefix, "_program_rank.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
utils::write.table(coef_rows, paste0(table_prefix, "_regulator_coefficients.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
utils::write.table(permutation_rows, paste0(table_prefix, "_permutation.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

if (plot_enabled) {
  ensure_parent_dir(paste0(plot_prefix, ".pdf"))
  plot_df <- graph_long
  if (nrow(plot_df) == 0) {
    grDevices::pdf(paste0(plot_prefix, ".pdf"), width = 8, height = 5)
    plot.new()
    title(main = sprintf("%s: no concordant graph genes", trait_id))
    grDevices::dev.off()
    grDevices::png(paste0(plot_prefix, ".png"), width = 8, height = 5, units = "in", res = 300)
    plot.new()
    title(main = sprintf("%s: no concordant graph genes", trait_id))
    grDevices::dev.off()
  } else {
    suppressPackageStartupMessages(library(ggplot2))
    plot_df$Program <- factor(plot_df$Program, levels = unique(selected_rows$Program))
    plot_df$side <- factor(plot_df$side, levels = c("program_loading", "regulator"))
    g <- ggplot(plot_df, aes(x = side, y = reorder(gene, abs_gamma), fill = post_mean))
    g <- g + geom_tile(color = "white", linewidth = 0.2)
    g <- g + facet_wrap(~ Program, scales = "free_y", ncol = 2)
    g <- g + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0)
    g <- g + theme_minimal(base_size = 12, base_family = "Helvetica")
    g <- g + theme(axis.title = element_blank(), panel.grid = element_blank())
    g <- g + labs(fill = "gamma", title = sprintf("%s Trait-Program-Gene graph genes", trait_id))
    grDevices::pdf(paste0(plot_prefix, ".pdf"), width = 10, height = max(6, 0.28 * nrow(plot_df) + 3))
    print(g)
    grDevices::dev.off()
    grDevices::png(paste0(plot_prefix, ".png"), width = 10, height = max(6, 0.28 * nrow(plot_df) + 3), units = "in", res = 300)
    print(g)
    grDevices::dev.off()
  }
}

quit(save = "no")
