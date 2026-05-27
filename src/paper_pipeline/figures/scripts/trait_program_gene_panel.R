args <- commandArgs(trailingOnly = TRUE)

program_association_dir <- as.character(args[1])
regulation_dir <- as.character(args[2])
spectra_path <- as.character(args[3])
gene_map_path <- as.character(args[4])
posterior_path <- as.character(args[5])
trait_id <- as.character(args[6])
k <- as.numeric(args[7])
table_prefix <- as.character(args[8])
plot_prefix <- as.character(args[9])
shet_path <- if (length(args) >= 10) as.character(args[10]) else ""
program_n <- if (length(args) >= 11) as.numeric(args[11]) else 5
regulator_n <- if (length(args) >= 12) as.numeric(args[12]) else 3
max_genes_per_side <- if (length(args) >= 13) as.numeric(args[13]) else 8
hit_abs_gamma_threshold <- if (length(args) >= 14) as.numeric(args[14]) else 0.1
loading_top_n <- if (length(args) >= 15) as.numeric(args[15]) else 200
regulator_fdr_threshold <- if (length(args) >= 16) as.numeric(args[16]) else 0.05
random_iterations <- if (length(args) >= 17) as.numeric(args[17]) else 10000
render_plot <- if (length(args) >= 18) as.character(args[18]) else "1"
association_trait_file <- if (length(args) >= 19) as.character(args[19]) else ""

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "helpers.R"))

options(scipen = 10)

plot_enabled <- render_plot %in% c("1", "true", "TRUE", "yes", "YES")

make_empty_side_hits <- function() {
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
    rank_within_side = integer(),
    program_trait_sign = character(),
    regulator_program_sign = character(),
    predicted_sign = character(),
    post_mean_sign = character(),
    is_concordant = logical(),
    is_discordant = logical(),
    display_bucket = character(),
    display_bucket_label = character(),
    display_column = character(),
    display_column_rank = integer(),
    program_label = character(),
    gene_label = character(),
    display_rank = integer(),
    x = numeric(),
    y = numeric(),
    panel_row = integer(),
    y_global = numeric(),
    has_overlap = logical(),
    empty_reason = character(),
    stringsAsFactors = FALSE
  )
}

sign_to_label <- function(values) {
  values <- as.numeric(values)
  labels <- rep("zero", length(values))
  labels[is.na(values)] <- ""
  labels[!is.na(values) & values > 0] <- "positive"
  labels[!is.na(values) & values < 0] <- "negative"
  labels
}

write_placeholder_plot <- function(program_summary, plot_prefix, trait_id, max_genes_per_side, message_text) {
  ensure_parent_dir(paste0(plot_prefix, ".pdf"))
  library(ggplot2)

  panel_rows <- nrow(program_summary)
  gene_text_size <- max(3.2, min(5.2, 5.8 - max_genes_per_side * 0.22))
  plot_height <- max(5, panel_rows * 1.1 + 3.5)
  title_y <- if (panel_rows > 0) max(program_summary$y_center) + 2 else 4
  y_max <- if (panel_rows > 0) title_y + 1.5 else 6

  g <- ggplot()
  if (panel_rows > 0) {
    g <- g + geom_segment(
      data = program_summary,
      aes(x = -0.42, xend = 0.42, y = y_center, yend = y_center),
      linewidth = 0.4,
      color = "grey88"
    )
    g <- g + geom_label(
      data = program_summary,
      aes(x = 0, y = y_center, label = program_label, fill = color),
      label.size = 0.15,
      label.padding = grid::unit(0.18, "lines"),
      size = gene_text_size * 0.72,
      fontface = "bold",
      color = "black"
    )
  }

  g <- g + annotate(
    "text",
    x = 0,
    y = title_y,
    label = trait_id,
    hjust = 0.5,
    fontface = "bold",
    size = gene_text_size * 1.05
  )
  g <- g + annotate(
    "text",
    x = 0,
    y = title_y - 1.05,
    label = message_text,
    hjust = 0.5,
    size = gene_text_size * 0.78,
    color = "grey35"
  )
  g <- g + scale_fill_manual(
    values = c(
      "other" = "#ECECEC",
      "program_enriched" = "#F9D48B",
      "regulator_enriched" = "#A9C7E6",
      "both_enriched" = "#A9D6A4"
    ),
    drop = FALSE
  )
  g <- g + coord_cartesian(xlim = c(-0.82, 0.82), ylim = c(0, y_max), clip = "off")
  g <- g + theme_void(base_family = "Helvetica")
  g <- g + theme(
    legend.position = "none",
    plot.margin = margin(18, 28, 18, 28)
  )

  ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = 10, height = plot_height)
  ggsave(plot = g, filename = paste0(plot_prefix, ".png"), width = 10, height = plot_height, dpi = 300)
}

program_association_files_exist <- function(association_dir, trait_file, k) {
  all(file.exists(c(
    file.path(association_dir, paste0("programs_enrichment_K", k, "_", trait_file)),
    file.path(association_dir, paste0("regulators_enrichment_K", k, "_", trait_file))
  )))
}

resolve_association_trait_file <- function(association_dir, association_trait_file, posterior_path, k) {
  posterior_file <- basename(posterior_path)
  posterior_as_csv <- sub("\\.per_gene_estimates\\.tsv$", ".csv", posterior_file)
  posterior_as_tsv <- sub("\\.per_gene_estimates\\.tsv$", ".tsv", posterior_file)
  candidates <- unique(c(association_trait_file, posterior_file, posterior_as_csv, posterior_as_tsv))
  candidates <- candidates[nzchar(candidates)]

  if (length(candidates) == 0) {
    return(posterior_file)
  }

  for (candidate in candidates) {
    if (program_association_files_exist(association_dir, candidate, k)) {
      if (!identical(candidate, association_trait_file)) {
        message(sprintf(
          "Using fallback ProgramLevel trait file '%s' because requested '%s' was not found",
          candidate,
          association_trait_file
        ))
      }
      return(candidate)
    }
  }

  association_trait_file
}

annotate_program_summary <- function(program_summary) {
  program_summary$Program <- normalize_program_ids(program_summary$Program)
  bonferroni_cutoff <- if (nrow(program_summary) > 0) 0.05 / nrow(program_summary) else Inf

  program_summary$program_sig <- !is.na(program_summary$MEANgamma_top100_shet_adjusted_P) &
    program_summary$MEANgamma_top100_shet_adjusted_P < bonferroni_cutoff
  program_summary$regulator_sig <- !is.na(program_summary$P_withShet) &
    program_summary$P_withShet < bonferroni_cutoff
  program_summary$priority_tier <- ifelse(
    program_summary$program_sig & program_summary$regulator_sig,
    1,
    ifelse(program_summary$program_sig | program_summary$regulator_sig, 2, 3)
  )
  program_summary$priority_score <- pmax(abs(program_summary$program_score), abs(program_summary$regulator_score), na.rm = TRUE)
  program_summary$priority_score[is.na(program_summary$priority_score)] <- 0
  program_summary$color <- ifelse(
    program_summary$program_sig & program_summary$regulator_sig,
    "both_enriched",
    ifelse(
      program_summary$program_sig,
      "program_enriched",
      ifelse(program_summary$regulator_sig, "regulator_enriched", "other")
    )
  )
  program_summary$color <- factor(
    program_summary$color,
    levels = c("other", "program_enriched", "regulator_enriched", "both_enriched")
  )
  program_summary
}

posterior_df <- data.table::fread(posterior_path, data.table = FALSE)
if (!all(c("ensg", "post_mean") %in% names(posterior_df))) {
  stop(sprintf("Posterior file %s must have columns: ensg, post_mean", posterior_path))
}
posterior_df$ensg <- normalize_ensg_ids(posterior_df$ensg)
posterior_df$post_mean[is.infinite(posterior_df$post_mean)] <- NA

lookups <- build_gene_id_lookups(gene_map_path)
posterior_df$gene <- unname(lookups$gene_lookup[posterior_df$ensg])
posterior_df$gene[is.na(posterior_df$gene) | !nzchar(posterior_df$gene)] <- posterior_df$ensg[is.na(posterior_df$gene) | !nzchar(posterior_df$gene)]
posterior_df$gamma_sign <- ifelse(posterior_df$post_mean > 0, "positive", ifelse(posterior_df$post_mean < 0, "negative", "zero"))
posterior_df$abs_gamma <- abs(posterior_df$post_mean)

trait_file <- resolve_association_trait_file(program_association_dir, association_trait_file, posterior_path, k)
message(sprintf("Using ProgramLevel association trait file: %s", trait_file))
message(sprintf("Using posterior file: %s", basename(posterior_path)))
program_summary_all <- read_program_regulator_summary(program_association_dir, trait_file, k, character())
program_summary_all <- annotate_program_summary(program_summary_all)

if (!nzchar(shet_path) || !file.exists(shet_path)) {
  stop(sprintf("trait_program_gene_panel.R requires shet_path, got: %s", shet_path))
}

program_n <- max(1, as.integer(program_n))
regulator_n <- max(1, as.integer(regulator_n))
random_iterations <- max(100, as.integer(random_iterations))

compute_program_selection <- function(posterior_df, spectra_path, gene_map_path, shet_path, k, top_n, random_iterations, n_select) {
  spectra_raw <- data.table::fread(spectra_path, data.table = FALSE)
  if (ncol(spectra_raw) < 1) {
    return(data.frame(Program = character(), P = numeric(), meanG = numeric(), stringsAsFactors = FALSE))
  }
  if (is.character(spectra_raw[[1]])) {
    row.names(spectra_raw) <- spectra_raw[[1]]
    spectra_raw <- spectra_raw[, -1, drop = FALSE]
  }
  spectra_mat <- t(as.matrix(spectra_raw))
  actual_k <- min(k, ncol(spectra_mat))
  colnames(spectra_mat) <- paste0("P", seq_len(ncol(spectra_mat)))

  lookups <- build_gene_id_lookups(gene_map_path)
  spectra_ensg <- normalize_ensg_ids(row.names(spectra_mat))
  spectra_genes <- unname(lookups$gene_lookup[spectra_ensg])
  missing_gene <- is.na(spectra_genes) | !nzchar(spectra_genes)
  spectra_genes[missing_gene] <- spectra_ensg[missing_gene]

  lof_df <- posterior_df[, c("ensg", "post_mean"), drop = FALSE]
  lof_df$ensg <- normalize_ensg_ids(lof_df$ensg)
  lof_df$gene <- unname(lookups$gene_lookup[lof_df$ensg])
  lof_df$gene[is.na(lof_df$gene) | !nzchar(lof_df$gene)] <- lof_df$ensg[is.na(lof_df$gene) | !nzchar(lof_df$gene)]
  lof_df <- lof_df[!is.na(lof_df$post_mean), , drop = FALSE]

  shet_df <- utils::read.table(shet_path, header = TRUE, stringsAsFactors = FALSE)
  shet_df$ensg <- normalize_ensg_ids(shet_df$ensg)
  shet_df <- shet_df[shet_df$ensg %in% lof_df$ensg, , drop = FALSE]

  rows <- list()
  for (program_index in seq_len(actual_k)) {
    program_name <- paste0("P", program_index)
    tmp <- data.frame(
      gene = spectra_genes,
      GEPscore = as.numeric(spectra_mat[, program_index]),
      ensg = spectra_ensg,
      stringsAsFactors = FALSE
    )
    tmp <- tmp[order(tmp$GEPscore, decreasing = TRUE, na.last = NA), , drop = FALSE]
    tmp <- tmp[!duplicated(tmp$ensg) & nzchar(tmp$ensg), , drop = FALSE]
    df <- merge(tmp, lof_df, by = "ensg")
    if (nrow(df) == 0) {
      next
    }
    df2 <- head(df[order(df$GEPscore, decreasing = TRUE), , drop = FALSE], top_n)
    if (nrow(df2) == 0) {
      next
    }
    mean_gamma_top <- mean(df2$post_mean)
    shet_tmp <- shet_df[shet_df$ensg %in% df2$ensg, , drop = FALSE]
    shet_tmp_all <- shet_df[shet_df$ensg %in% df$ensg, , drop = FALSE]
    if (nrow(shet_tmp) == 0 || nrow(shet_tmp_all) == 0 || !"shet_BIN" %in% names(shet_tmp_all)) {
      next
    }
    A <- table(shet_tmp$shet_BIN)
    random_means <- numeric()
    for (seed in seq_len(random_iterations)) {
      set.seed(seed)
      sampled_genes <- character()
      valid <- TRUE
      for (bin_name in names(A)) {
        candidates <- shet_tmp_all$ensg[shet_tmp_all$shet_BIN %in% bin_name]
        if (length(candidates) < A[[bin_name]]) {
          valid <- FALSE
          break
        }
        sampled_genes <- c(sampled_genes, sample(candidates, A[[bin_name]]))
      }
      if (!valid) {
        next
      }
      random_means <- c(random_means, mean(df$post_mean[df$ensg %in% sampled_genes], na.rm = TRUE))
    }
    if (length(random_means) == 0) {
      next
    }
    random_mean <- mean(random_means)
    rank_value <- rank(c(mean_gamma_top, random_means))[1]
    p1 <- (rank_value / length(random_means)) * 2
    p2 <- ((length(random_means) + 2 - rank_value) / length(random_means)) * 2
    rows[[length(rows) + 1]] <- data.frame(
      Program = program_name,
      P = min(p1, p2),
      meanG = mean_gamma_top - random_mean,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) {
    return(data.frame(Program = character(), P = numeric(), meanG = numeric(), stringsAsFactors = FALSE))
  }
  program_p_sum <- do.call(rbind, rows)
  program_p_sum <- program_p_sum[order(abs(program_p_sum$meanG), decreasing = TRUE), , drop = FALSE]
  program_p_sum <- program_p_sum[order(program_p_sum$P), , drop = FALSE]
  head(program_p_sum, n_select)
}

compute_regulator_selection <- function(posterior_df, regulation_dir, gene_map_path, shet_path, k, n_select) {
  lookups <- build_gene_id_lookups(gene_map_path)
  reg_mats <- read_regulation_matrices(regulation_dir, k)
  beta_df <- reg_mats$beta
  if (nrow(beta_df) == 0) {
    return(character())
  }

  lof_df <- posterior_df[, c("ensg", "post_mean"), drop = FALSE]
  lof_df$ensg <- normalize_ensg_ids(lof_df$ensg)
  lof_df$GENE <- unname(lookups$gene_lookup[lof_df$ensg])
  lof_df$GENE[is.na(lof_df$GENE) | !nzchar(lof_df$GENE)] <- lof_df$ensg[is.na(lof_df$GENE) | !nzchar(lof_df$GENE)]
  lof_df <- lof_df[!is.na(lof_df$post_mean), c("GENE", "ensg", "post_mean"), drop = FALSE]

  shet_df <- utils::read.table(shet_path, header = TRUE, stringsAsFactors = FALSE)
  shet_df$ensg <- normalize_ensg_ids(shet_df$ensg)
  shet_sub <- shet_df[shet_df$ensg %in% lof_df$ensg, , drop = FALSE]

  df <- merge(lof_df, beta_df, by = "GENE")
  df <- merge(df, shet_sub[, c("ensg", "shet"), drop = FALSE], by = "ensg")
  if (nrow(df) == 0) {
    return(character())
  }
  keep_cols <- c("post_mean", paste0("P", seq_len(k)), "shet")
  df <- df[, keep_cols, drop = FALSE]
  for (i in seq_len(ncol(df))) {
    finite_values <- df[[i]][is.finite(df[[i]])]
    if (length(finite_values) == 0) {
      next
    }
    df[[i]][is.infinite(df[[i]])] <- max(finite_values)
  }

  suppressPackageStartupMessages(library(leaps))
  fit <- tryCatch(
    leaps::regsubsets(post_mean ~ ., data = df, nbest = 1, nvmax = n_select + 1, really.big = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(character())
  }
  fit_sum <- summary(fit)[[1]]
  regulator_selected <- colnames(fit_sum)[fit_sum[nrow(fit_sum), ]]
  regulator_selected <- regulator_selected[!regulator_selected %in% c("(Intercept)", "shet")]
  if (length(regulator_selected) != n_select) {
    fit <- tryCatch(
      leaps::regsubsets(post_mean ~ ., data = df, nbest = 1, nvmax = n_select, really.big = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(normalize_program_ids(regulator_selected))
    }
    fit_sum <- summary(fit)[[1]]
    regulator_selected <- colnames(fit_sum)[fit_sum[nrow(fit_sum), ]]
    regulator_selected <- regulator_selected[!regulator_selected %in% c("(Intercept)", "shet")]
  }
  normalize_program_ids(regulator_selected)
}

program_selection_df <- compute_program_selection(
  posterior_df = posterior_df,
  spectra_path = spectra_path,
  gene_map_path = gene_map_path,
  shet_path = shet_path,
  k = k,
  top_n = loading_top_n,
  random_iterations = random_iterations,
  n_select = program_n
)
program_selected <- normalize_program_ids(program_selection_df$Program)
regulator_selected <- compute_regulator_selection(
  posterior_df = posterior_df,
  regulation_dir = regulation_dir,
  gene_map_path = gene_map_path,
  shet_path = shet_path,
  k = k,
  n_select = regulator_n
)
selected_programs <- unique(c(program_selected, regulator_selected))

program_summary <- program_summary_all[
  match(selected_programs, program_summary_all$Program, nomatch = 0),
  ,
  drop = FALSE
]
program_summary$selected_by_program <- program_summary$Program %in% program_selected
program_summary$selected_by_regulator <- program_summary$Program %in% regulator_selected
program_summary$program_selection_p <- unname(stats::setNames(program_selection_df$P, program_selection_df$Program)[program_summary$Program])
program_summary$program_selection_meanG <- unname(stats::setNames(program_selection_df$meanG, program_selection_df$Program)[program_summary$Program])

spectra_top <- read_spectra_top_genes(spectra_path, gene_map_path, k, loading_top_n)
spectra_top <- spectra_top[spectra_top$Program %in% selected_programs, , drop = FALSE]

regulator_df <- read_regulator_long(regulation_dir, k)
regulator_df <- regulator_df[regulator_df$Program %in% selected_programs, , drop = FALSE]
regulator_df <- regulator_df[regulator_df$fdr <= regulator_fdr_threshold, , drop = FALSE]
regulator_df <- regulator_df[order(regulator_df$Program, regulator_df$fdr, -abs(regulator_df$beta), regulator_df$gene), , drop = FALSE]
if (nrow(regulator_df) > 0) {
  regulator_df$gene_ensg <- lookups$resolve_to_ensg(regulator_df$gene)
  regulator_df$gene_symbol <- unname(lookups$gene_lookup[regulator_df$gene_ensg])
  regulator_df$gene_symbol[is.na(regulator_df$gene_symbol) | !nzchar(regulator_df$gene_symbol)] <- regulator_df$gene[
    is.na(regulator_df$gene_symbol) | !nzchar(regulator_df$gene_symbol)
  ]
  regulator_df$gene <- regulator_df$gene_symbol
  regulator_df$rank_within_side <- ave(
    regulator_df$fdr,
    regulator_df$Program,
    FUN = function(x) seq_along(x)
  )
} else {
  regulator_df$rank_within_side <- numeric()
}

trait_hits <- posterior_df[!is.na(posterior_df$abs_gamma) & posterior_df$abs_gamma > hit_abs_gamma_threshold, , drop = FALSE]

side_hit_cols <- c("Program", "side", "gene", "ensg", "post_mean", "abs_gamma", "gamma_sign", "membership_score", "rank_within_side")

loading_hits <- merge(
  spectra_top,
  trait_hits[, c("ensg", "gene", "post_mean", "abs_gamma", "gamma_sign")],
  by = "gene",
  all = FALSE
)
if (nrow(loading_hits) > 0) {
  loading_hits$side <- "program_loading"
  loading_hits$membership_score <- loading_hits$weight
  loading_hits$effect_score <- loading_hits$post_mean
  loading_hits$effect_rank <- ave(
    -loading_hits$abs_gamma,
    loading_hits$Program,
    FUN = rank,
    ties.method = "first"
  )
  loading_hits <- loading_hits[order(loading_hits$Program, loading_hits$effect_rank, loading_hits$rank_within_side), , drop = FALSE]
  loading_hits <- do.call(
    rbind,
    lapply(split(loading_hits, loading_hits$Program), function(df) head(df, max_genes_per_side))
  )
} else {
  loading_hits <- make_empty_side_hits()[, side_hit_cols, drop = FALSE]
}

regulator_hits <- merge(
  regulator_df,
  trait_hits[, c("ensg", "gene", "post_mean", "abs_gamma", "gamma_sign")],
  by = "gene",
  all = FALSE
)
if (nrow(regulator_hits) > 0) {
  regulator_hits$side <- "regulator"
  regulator_hits$membership_score <- regulator_hits$beta
  regulator_hits$effect_score <- regulator_hits$post_mean
  regulator_hits$effect_rank <- ave(
    -regulator_hits$abs_gamma,
    regulator_hits$Program,
    FUN = rank,
    ties.method = "first"
  )
  regulator_hits <- regulator_hits[order(regulator_hits$Program, regulator_hits$effect_rank, regulator_hits$rank_within_side), , drop = FALSE]
  regulator_hits <- do.call(
    rbind,
    lapply(split(regulator_hits, regulator_hits$Program), function(df) head(df, max_genes_per_side))
  )
} else {
  regulator_hits <- make_empty_side_hits()[, side_hit_cols, drop = FALSE]
}

side_hit_parts <- list()
if (nrow(loading_hits) > 0) {
  side_hit_parts[[length(side_hit_parts) + 1]] <- loading_hits[, side_hit_cols]
}
if (nrow(regulator_hits) > 0) {
  side_hit_parts[[length(side_hit_parts) + 1]] <- regulator_hits[, side_hit_cols]
}
if (length(side_hit_parts) > 0) {
  side_hits <- do.call(rbind, side_hit_parts)
} else {
  side_hits <- make_empty_side_hits()[, side_hit_cols, drop = FALSE]
}

has_overlap <- nrow(side_hits) > 0
empty_reason <- ""
if (!has_overlap) {
  empty_reason <- if (length(selected_programs) == 0) "no_selected_programs" else "no_trait_program_gene_overlap"
  message(sprintf("No trait-program gene overlaps for %s; writing empty outputs (%s)", trait_id, empty_reason))
}

program_summary$Program <- as.character(program_summary$Program)
if (nrow(program_summary) > 0 && length(selected_programs) > 0) {
  program_summary <- program_summary[match(selected_programs, program_summary$Program, nomatch = 0), , drop = FALSE]
}

loading_gene_count <- integer(nrow(program_summary))
regulator_gene_count <- integer(nrow(program_summary))
if (has_overlap && nrow(program_summary) > 0) {
  loading_lookup <- table(side_hits$Program[side_hits$side == "program_loading"])
  regulator_lookup <- table(side_hits$Program[side_hits$side == "regulator"])
  loading_gene_count <- as.integer(loading_lookup[program_summary$Program])
  regulator_gene_count <- as.integer(regulator_lookup[program_summary$Program])
  loading_gene_count[is.na(loading_gene_count)] <- 0L
  regulator_gene_count[is.na(regulator_gene_count)] <- 0L
}
program_summary$loading_gene_count <- loading_gene_count
program_summary$regulator_gene_count <- regulator_gene_count
program_summary$program_label <- sprintf(
  "%s  L:%d  R:%d",
  program_summary$Program,
  program_summary$loading_gene_count,
  program_summary$regulator_gene_count
)
if (!"selected_by_program" %in% names(program_summary)) {
  program_summary$selected_by_program <- program_summary$Program %in% program_selected
}
if (!"selected_by_regulator" %in% names(program_summary)) {
  program_summary$selected_by_regulator <- program_summary$Program %in% regulator_selected
}
program_summary$trait_id <- trait_id
program_summary$has_overlap <- has_overlap
program_summary$empty_reason <- if (has_overlap) "" else empty_reason

panel_rows <- nrow(program_summary)
row_spacing <- max_genes_per_side + 3
program_summary$panel_row <- if (panel_rows > 0) seq_len(panel_rows) else integer()
program_summary$y_center <- if (panel_rows > 0) {
  (panel_rows - program_summary$panel_row) * row_spacing + (max_genes_per_side / 2)
} else {
  numeric()
}
program_mean_sign <- sign(program_summary$MEANgamma_top100 - program_summary$shet_adjusted_random_mean)
names(program_mean_sign) <- program_summary$Program
program_summary$program_trait_sign <- sign_to_label(program_mean_sign)

if (has_overlap) {
  regulator_sign_lookup <- if (nrow(regulator_df) > 0) {
    setNames(sign(regulator_df$beta), paste(regulator_df$Program, regulator_df$gene, sep = "::"))
  } else {
    character()
  }

  is_loading <- side_hits$side == "program_loading"
  side_hits$program_trait_sign_value <- unname(program_mean_sign[side_hits$Program])
  side_hits$program_trait_sign_value[is.na(side_hits$program_trait_sign_value)] <- 0
  side_hits$regulator_program_sign_value <- 0
  if (any(!is_loading)) {
    reg_keys <- paste(side_hits$Program[!is_loading], side_hits$gene[!is_loading], sep = "::")
    side_hits$regulator_program_sign_value[!is_loading] <- unname(regulator_sign_lookup[reg_keys])
  }
  side_hits$regulator_program_sign_value[is.na(side_hits$regulator_program_sign_value)] <- 0
  side_hits$predicted_sign_value <- side_hits$program_trait_sign_value
  side_hits$predicted_sign_value[!is_loading] <- (
    side_hits$program_trait_sign_value[!is_loading] *
      side_hits$regulator_program_sign_value[!is_loading]
  )
  side_hits$post_mean_sign_value <- sign(side_hits$post_mean)
  side_hits$is_concordant <- (
    side_hits$predicted_sign_value != 0 &
      side_hits$post_mean_sign_value != 0 &
      side_hits$predicted_sign_value == side_hits$post_mean_sign_value
  )
  side_hits$is_discordant <- (
    side_hits$predicted_sign_value != 0 &
      side_hits$post_mean_sign_value != 0 &
      side_hits$predicted_sign_value == -side_hits$post_mean_sign_value
  )
  side_hits$program_trait_sign <- sign_to_label(side_hits$program_trait_sign_value)
  side_hits$regulator_program_sign <- ""
  side_hits$regulator_program_sign[!is_loading] <- sign_to_label(side_hits$regulator_program_sign_value[!is_loading])
  side_hits$predicted_sign <- sign_to_label(side_hits$predicted_sign_value)
  side_hits$post_mean_sign <- sign_to_label(side_hits$post_mean_sign_value)
  side_hits$display_bucket <- ifelse(
    is_loading,
    "program_genes",
    ifelse(side_hits$regulator_program_sign_value < 0, "negative_regulators", "positive_regulators")
  )
  side_hits$display_bucket_label <- ifelse(
    side_hits$display_bucket == "program_genes",
    "Program genes",
    ifelse(side_hits$display_bucket == "negative_regulators", "Negative regulators", "Positive regulators")
  )
  side_hits$display_column <- ifelse(
    side_hits$predicted_sign_value == 0,
    ifelse(side_hits$post_mean_sign_value < 0, "right", "left"),
    ifelse(side_hits$post_mean_sign_value == side_hits$predicted_sign_value, "left", "right")
  )

  label_lookup <- setNames(program_summary$program_label, program_summary$Program)
  side_hits$program_label <- unname(label_lookup[side_hits$Program])

  side_hits$gene_label <- ifelse(side_hits$is_discordant, paste0("(", side_hits$gene, ")"), side_hits$gene)
  max_label_chars <- max(12, min(22, floor(90 / max(1, max_genes_per_side))))
  side_hits$gene_label <- truncate_gene_label(side_hits$gene_label, max_chars = max_label_chars)

  side_hits$display_rank <- ave(
    side_hits$rank_within_side,
    side_hits$Program,
    side_hits$side,
    FUN = function(x) seq_along(x)
  )
  side_hits$display_column_rank <- ave(
    side_hits$display_rank,
    side_hits$Program,
    side_hits$display_bucket,
    side_hits$display_column,
    FUN = function(x) seq_along(x)
  )
  side_hits$x <- ifelse(side_hits$side == "program_loading", -1, 1)
  side_hits$y <- NA_real_
  program_groups <- split(seq_len(nrow(side_hits)), side_hits$Program, drop = TRUE)
  for (idx in program_groups) {
    df <- side_hits[idx, , drop = FALSE]
    loading_idx <- idx[df$side == "program_loading"]
    regulator_idx <- idx[df$side == "regulator"]
    if (length(loading_idx) > 0) {
      side_hits$y[loading_idx] <- rev(seq_len(length(loading_idx)))
    }
    if (length(regulator_idx) > 0) {
      side_hits$y[regulator_idx] <- rev(seq_len(length(regulator_idx)))
    }
  }

  program_row_lookup <- setNames(seq_len(panel_rows), program_summary$Program)
  side_hits$panel_row <- unname(program_row_lookup[side_hits$Program])
  side_hits$y_global <- (panel_rows - side_hits$panel_row) * row_spacing + side_hits$y
  side_hits$trait_id <- trait_id
  side_hits$has_overlap <- TRUE
  side_hits$empty_reason <- ""
  side_hits <- side_hits[, names(make_empty_side_hits()), drop = FALSE]
  side_hits$membership_score <- as.numeric(side_hits$membership_score)
} else {
  side_hits <- make_empty_side_hits()
}

table_long_path <- paste0(table_prefix, "_long.tsv")
table_program_path <- paste0(table_prefix, "_programs.tsv")
ensure_parent_dir(table_long_path)
utils::write.table(side_hits, table_long_path, row.names = FALSE, sep = "\t", quote = FALSE)
utils::write.table(program_summary, table_program_path, row.names = FALSE, sep = "\t", quote = FALSE)

if (!plot_enabled) {
  quit(save = "no")
}

if (!has_overlap) {
  placeholder_message <- if (identical(empty_reason, "no_selected_programs")) {
    "No programs passed current selection filters"
  } else {
    "No overlapping genes passed current filters"
  }
  write_placeholder_plot(program_summary, plot_prefix, trait_id, max_genes_per_side, placeholder_message)
  quit(save = "no")
}

ensure_parent_dir(paste0(plot_prefix, ".pdf"))
library(ggplot2)
library(grid)

draw_model_panel <- function(
    program_summary,
    side_hits,
    trait_id,
    plot_prefix,
    max_genes_per_side,
    hit_abs_gamma_threshold,
    loading_top_n,
    regulator_fdr_threshold
) {
  clean_label <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("[\r\n\t]+", " ", x)
    x <- gsub("\\s+", " ", x)
    trimws(x)
  }

  format_signed_score <- function(values) {
    values <- as.numeric(values)
    out <- ifelse(
      is.na(values),
      "NA",
      ifelse(
        is.finite(values),
        sprintf("%+.1f", values),
        sprintf("%sInf", ifelse(values > 0, "+", "-"))
      )
    )
    out
  }

  signed_value_color <- function(values, max_abs) {
    values <- as.numeric(values)
    values[!is.finite(values)] <- 0
    max_abs <- max(as.numeric(max_abs), 1e-12)

    pos_palette <- grDevices::colorRampPalette(c("#F2B6A7", "#B2182B"))(101)
    neg_palette <- grDevices::colorRampPalette(c("#B8D7EA", "#2166AC"))(101)
    intensity <- pmin(1, abs(values) / max_abs)
    idx <- pmax(1, pmin(101, floor(intensity * 100) + 1))

    colors <- rep("#7A808A", length(values))
    colors[values > 0] <- pos_palette[idx[values > 0]]
    colors[values < 0] <- neg_palette[idx[values < 0]]
    colors
  }

  split_gene_lines <- function(df, n = max_genes_per_side) {
    if (nrow(df) == 0) {
      return(data.frame(
        gene_label = character(),
        post_mean = numeric(),
        is_discordant = logical(),
        stringsAsFactors = FALSE
      ))
    }
    df <- df[order(df$display_rank, -df$abs_gamma, df$gene), , drop = FALSE]
    df <- head(df, n)
    data.frame(
      gene_label = clean_label(df$gene_label),
      post_mean = as.numeric(df$post_mean),
      is_discordant = as.logical(df$is_discordant),
      stringsAsFactors = FALSE
    )
  }

  draw_gene_list <- function(lines, x, y_top, line_gap, hjust, font_size, gamma_max) {
    if (nrow(lines) == 0) {
      grid.text(
        "no overlap",
        x = unit(x, "npc"),
        y = unit(y_top, "npc"),
        just = c(ifelse(hjust == 0, "left", "right"), "top"),
        gp = gpar(col = "#9CA3AF", fontsize = font_size * 0.82, fontface = "italic")
      )
      return(invisible(NULL))
    }

    colors <- signed_value_color(lines$post_mean, gamma_max)
    for (i in seq_len(nrow(lines))) {
      label <- clean_label(lines$gene_label[i])
      fontface <- if (isTRUE(lines$is_discordant[i])) "bold.italic" else "plain"
      grid.text(
        label,
        x = unit(x, "npc"),
        y = unit(y_top - (i - 1) * line_gap, "npc"),
        just = c(ifelse(hjust == 0, "left", "right"), "top"),
        gp = gpar(col = colors[i], fontsize = font_size, fontface = fontface)
      )
    }
  }

  draw_arrow <- function(x0, x1, y, direction, col = "#4B5563", lwd = 2.1) {
    direction <- sign(as.numeric(direction))
    if (!is.finite(direction)) {
      direction <- 0
    }
    if (direction >= 0) {
      grid.segments(
        x0 = unit(x0, "npc"),
        x1 = unit(x1, "npc"),
        y0 = unit(y, "npc"),
        y1 = unit(y, "npc"),
        gp = gpar(col = col, lwd = lwd, lineend = "round"),
        arrow = arrow(length = unit(0.09, "inches"), type = "closed")
      )
    } else {
      grid.segments(
        x0 = unit(x0, "npc"),
        x1 = unit(x1, "npc"),
        y0 = unit(y, "npc"),
        y1 = unit(y, "npc"),
        gp = gpar(col = col, lwd = lwd, lineend = "round")
      )
      grid.segments(
        x0 = unit(x1, "npc"),
        x1 = unit(x1, "npc"),
        y0 = unit(y - 0.010, "npc"),
        y1 = unit(y + 0.010, "npc"),
        gp = gpar(col = col, lwd = lwd, lineend = "round")
      )
    }
  }

  draw_program_box <- function(x, y, label, score_label, fill, border = "#111827") {
    grid.roundrect(
      x = unit(x, "npc"),
      y = unit(y, "npc"),
      width = unit(0.060, "npc"),
      height = unit(0.032, "npc"),
      r = unit(0.006, "npc"),
      gp = gpar(fill = fill, col = border, lwd = 0.8)
    )
    grid.text(
      label,
      x = unit(x, "npc"),
      y = unit(y + 0.002, "npc"),
      gp = gpar(col = "#111827", fontsize = 9.2, fontface = "bold")
    )
    grid.text(
      score_label,
      x = unit(x, "npc"),
      y = unit(y - 0.031, "npc"),
      gp = gpar(col = "#4B5563", fontsize = 6.8, fontface = "bold")
    )
  }

  draw_gene_module <- function(x, y, width, height, title, lines, gamma_max, align = "left", fill = "#F8FAFC") {
    grid.roundrect(
      x = unit(x, "npc"),
      y = unit(y, "npc"),
      width = unit(width, "npc"),
      height = unit(height, "npc"),
      r = unit(0.008, "npc"),
      gp = gpar(fill = fill, col = "#E5E7EB", lwd = 0.6)
    )
    grid.text(
      title,
      x = unit(x - width / 2 + 0.014, "npc"),
      y = unit(y + height / 2 - 0.020, "npc"),
      just = c("left", "top"),
      gp = gpar(col = "#111827", fontsize = 7.4, fontface = "bold")
    )
    text_x <- if (identical(align, "right")) x + width / 2 - 0.014 else x - width / 2 + 0.014
    hjust <- if (identical(align, "right")) 1 else 0
    line_gap <- min(0.019, (height - 0.052) / max(1, max(nrow(lines), 1)))
    font_size <- max(5.1, min(7.2, 7.4 - max(0, nrow(lines) - 8) * 0.12))
    draw_gene_list(
      lines = lines,
      x = text_x,
      y_top = y + height / 2 - 0.047,
      line_gap = line_gap,
      hjust = hjust,
      font_size = font_size,
      gamma_max = gamma_max
    )
  }

  evidence_fill <- c(
    "other" = "#F3F4F6",
    "program_enriched" = "#F9D48B",
    "regulator_enriched" = "#A9C7E6",
    "both_enriched" = "#A9D6A4"
  )

  program_summary$Program <- as.character(program_summary$Program)
  side_hits$Program <- as.character(side_hits$Program)
  panel_programs <- program_summary[program_summary$selected_by_program, , drop = FALSE]
  panel_regulators <- program_summary[program_summary$selected_by_regulator, , drop = FALSE]
  panel_programs <- panel_programs[order(-abs(panel_programs$program_score), panel_programs$Program), , drop = FALSE]
  panel_regulators <- panel_regulators[order(-abs(panel_regulators$regulator_score), panel_regulators$Program), , drop = FALSE]
  gamma_max <- max(abs(as.numeric(side_hits$post_mean)), hit_abs_gamma_threshold, na.rm = TRUE)
  if (!is.finite(gamma_max)) {
    gamma_max <- max(hit_abs_gamma_threshold, 1)
  }

  n_left <- nrow(panel_programs)
  n_right <- nrow(panel_regulators)
  n_rows <- max(n_left, n_right, 1)
  canvas_width <- 15.8
  canvas_height <- max(7.0, min(14.0, 2.6 + n_rows * 1.35))
  ensure_parent_dir(paste0(plot_prefix, ".pdf"))

  draw_page <- function() {
    grid.newpage()
    grid.rect(gp = gpar(fill = "white", col = NA))

    grid.text(
      "Programs selected for modeling trait associations",
      x = unit(0.5, "npc"),
      y = unit(0.965, "npc"),
      gp = gpar(col = "#111827", fontsize = 18, fontface = "bold")
    )
    grid.text(
      sprintf(
        "Trait: %s   |   %d program-burden modules, %d regulator-burden modules. Gene color encodes signed gamma; parentheses mark discordant direction.",
        trait_id,
        n_left,
        n_right
      ),
      x = unit(0.5, "npc"),
      y = unit(0.936, "npc"),
      gp = gpar(col = "#4B5563", fontsize = 8.8)
    )

    grid.roundrect(
      x = unit(0.5, "npc"),
      y = unit(0.864, "npc"),
      width = unit(0.135, "npc"),
      height = unit(0.038, "npc"),
      r = unit(0.008, "npc"),
      gp = gpar(fill = "#111827", col = NA)
    )
    grid.text(
      trait_id,
      x = unit(0.5, "npc"),
      y = unit(0.866, "npc"),
      gp = gpar(col = "white", fontsize = 10.5, fontface = "bold")
    )

    grid.text(
      sprintf("Program burden selected\nprogram genes (top %d overlap)", loading_top_n),
      x = unit(0.247, "npc"),
      y = unit(0.806, "npc"),
      gp = gpar(col = "#111827", fontsize = 9.2, fontface = "bold", lineheight = 0.95)
    )
    grid.text(
      sprintf("Regulator-burden selected\nregulators (FDR <= %.3g)", regulator_fdr_threshold),
      x = unit(0.753, "npc"),
      y = unit(0.806, "npc"),
      gp = gpar(col = "#111827", fontsize = 9.2, fontface = "bold", lineheight = 0.95)
    )

    row_top <- 0.745
    row_bottom <- 0.155
    row_step <- if (n_rows <= 1) 0 else (row_top - row_bottom) / (n_rows - 1)
    module_h <- min(0.118, max(0.080, row_step * 0.73))

    if (n_left > 0) {
      for (i in seq_len(n_left)) {
        row <- panel_programs[i, , drop = FALSE]
        y <- if (n_rows <= 1) (row_top + row_bottom) / 2 else row_top - (i - 1) * row_step
        if (i %% 2 == 0) {
          grid.rect(
            x = unit(0.250, "npc"),
            y = unit(y, "npc"),
            width = unit(0.440, "npc"),
            height = unit(module_h * 1.16, "npc"),
            gp = gpar(fill = "#F9FAFB", col = NA)
          )
        }
        hits <- side_hits[side_hits$Program == row$Program & side_hits$side == "program_loading", , drop = FALSE]
        lines <- split_gene_lines(hits)
        draw_gene_module(
          x = 0.244,
          y = y,
          width = 0.270,
          height = module_h,
          title = sprintf("%s genes", row$Program),
          lines = lines,
          gamma_max = gamma_max,
          align = "right",
          fill = "#FFFFFF"
        )
        draw_arrow(
          x0 = 0.382,
          x1 = 0.462,
          y = y,
          direction = row$program_score,
          col = "#4B5563",
          lwd = 1.7 + 1.0 * min(1, abs(as.numeric(row$program_score)) / 4)
        )
        draw_program_box(
          x = 0.415,
          y = y + 0.030,
          label = row$Program,
          score_label = sprintf("B %s", format_signed_score(row$program_score)),
          fill = evidence_fill[as.character(row$color)]
        )
      }
    }

    if (n_right > 0) {
      for (i in seq_len(n_right)) {
        row <- panel_regulators[i, , drop = FALSE]
        y <- if (n_rows <= 1) (row_top + row_bottom) / 2 else row_top - (i - 1) * row_step
        if (i %% 2 == 0) {
          grid.rect(
            x = unit(0.750, "npc"),
            y = unit(y, "npc"),
            width = unit(0.440, "npc"),
            height = unit(module_h * 1.16, "npc"),
            gp = gpar(fill = "#F9FAFB", col = NA)
          )
        }
        hits <- side_hits[side_hits$Program == row$Program & side_hits$side == "regulator", , drop = FALSE]
        pos_hits <- hits[as.numeric(hits$membership_score) >= 0, , drop = FALSE]
        neg_hits <- hits[as.numeric(hits$membership_score) < 0, , drop = FALSE]
        pos_lines <- split_gene_lines(pos_hits, max(1, floor(max_genes_per_side / 2)))
        neg_lines <- split_gene_lines(neg_hits, max(1, ceiling(max_genes_per_side / 2)))
        if (nrow(pos_lines) == 0 && nrow(neg_lines) == 0) {
          all_lines <- split_gene_lines(hits)
          pos_lines <- all_lines
        }

        draw_program_box(
          x = 0.585,
          y = y + 0.030,
          label = row$Program,
          score_label = sprintf("R %s", format_signed_score(row$regulator_score)),
          fill = evidence_fill[as.character(row$color)]
        )
        draw_arrow(
          x0 = 0.538,
          x1 = 0.620,
          y = y,
          direction = row$regulator_score,
          col = "#4B5563",
          lwd = 1.7 + 1.0 * min(1, abs(as.numeric(row$regulator_score)) / 4)
        )
        draw_arrow(
          x0 = 0.620,
          x1 = 0.502,
          y = y,
          direction = row$regulator_score,
          col = "#9CA3AF",
          lwd = 0.9
        )

        draw_gene_module(
          x = 0.724,
          y = y + module_h * 0.24,
          width = 0.245,
          height = module_h * 0.48,
          title = "positive regulators",
          lines = pos_lines,
          gamma_max = gamma_max,
          align = "left",
          fill = "#FFF7ED"
        )
        draw_gene_module(
          x = 0.724,
          y = y - module_h * 0.30,
          width = 0.245,
          height = module_h * 0.48,
          title = "negative regulators",
          lines = neg_lines,
          gamma_max = gamma_max,
          align = "left",
          fill = "#EFF6FF"
        )
      }
    }

    legend_y <- 0.080
    grid.points(
      x = unit(c(0.160, 0.260, 0.360), "npc"),
      y = unit(rep(legend_y, 3), "npc"),
      pch = 16,
      size = unit(0.018, "npc"),
      gp = gpar(col = c("#2166AC", "#7A808A", "#B2182B"))
    )
    grid.text("negative gamma", x = unit(0.180, "npc"), y = unit(legend_y, "npc"), just = "left", gp = gpar(col = "#2166AC", fontsize = 7.4, fontface = "bold"))
    grid.text("near 0", x = unit(0.280, "npc"), y = unit(legend_y, "npc"), just = "left", gp = gpar(col = "#7A808A", fontsize = 7.4, fontface = "bold"))
    grid.text("positive gamma", x = unit(0.380, "npc"), y = unit(legend_y, "npc"), just = "left", gp = gpar(col = "#B2182B", fontsize = 7.4, fontface = "bold"))
    grid.text(
      "Arrow head = positive program/regulator direction; flat cap = negative direction.",
      x = unit(0.575, "npc"),
      y = unit(legend_y, "npc"),
      just = "left",
      gp = gpar(col = "#374151", fontsize = 7.4)
    )

    grid.text(
      sprintf(
        "Filters: |gamma| >= %.3g; program genes from top %d loading genes; regulators at FDR <= %.3g. B = program burden score; R = regulator-burden correlation score.",
        hit_abs_gamma_threshold,
        loading_top_n,
        regulator_fdr_threshold
      ),
      x = unit(0.020, "npc"),
      y = unit(0.032, "npc"),
      just = "left",
      gp = gpar(col = "#4B5563", fontsize = 7.2)
    )
  }

  draw_centered_page <- function() {
    grid.newpage()
    grid.rect(gp = gpar(fill = "white", col = NA))

    trait_x <- 0.500
    trait_y <- 0.500
    trait_w <- 0.150
    trait_h <- 0.062
    trait_left <- trait_x - trait_w / 2
    trait_right <- trait_x + trait_w / 2

    draw_trait_arrow <- function(x0, y0, x1, y1, direction, col = "#4B5563", lwd = 2.0) {
      direction <- sign(as.numeric(direction))
      if (!is.finite(direction)) {
        direction <- 0
      }
      if (direction >= 0) {
        grid.segments(
          x0 = unit(x0, "npc"),
          x1 = unit(x1, "npc"),
          y0 = unit(y0, "npc"),
          y1 = unit(y1, "npc"),
          gp = gpar(col = col, lwd = lwd, lineend = "round"),
          arrow = arrow(length = unit(0.085, "inches"), type = "closed")
        )
      } else {
        grid.segments(
          x0 = unit(x0, "npc"),
          x1 = unit(x1, "npc"),
          y0 = unit(y0, "npc"),
          y1 = unit(y1, "npc"),
          gp = gpar(col = col, lwd = lwd, lineend = "round")
        )
        grid.segments(
          x0 = unit(x1, "npc"),
          x1 = unit(x1, "npc"),
          y0 = unit(y1 - 0.012, "npc"),
          y1 = unit(y1 + 0.012, "npc"),
          gp = gpar(col = col, lwd = lwd, lineend = "round")
        )
      }
    }

    grid.text(
      "Programs selected for modeling trait associations",
      x = unit(0.5, "npc"),
      y = unit(0.965, "npc"),
      gp = gpar(col = "#111827", fontsize = 18, fontface = "bold")
    )
    grid.text(
      sprintf(
        "%d program-burden modules on the left; %d regulator-burden modules on the right. Gene color encodes signed gamma; parentheses mark discordant direction.",
        n_left,
        n_right
      ),
      x = unit(0.5, "npc"),
      y = unit(0.936, "npc"),
      gp = gpar(col = "#4B5563", fontsize = 8.8)
    )

    grid.text(
      sprintf("Program burden selected\nprogram loading genes (top %d overlap)", loading_top_n),
      x = unit(0.205, "npc"),
      y = unit(0.866, "npc"),
      gp = gpar(col = "#111827", fontsize = 9.2, fontface = "bold", lineheight = 0.95)
    )
    grid.text(
      sprintf("Regulator-burden selected\nregulator genes (FDR <= %.3g)", regulator_fdr_threshold),
      x = unit(0.795, "npc"),
      y = unit(0.866, "npc"),
      gp = gpar(col = "#111827", fontsize = 9.2, fontface = "bold", lineheight = 0.95)
    )

    grid.roundrect(
      x = unit(trait_x, "npc"),
      y = unit(trait_y, "npc"),
      width = unit(trait_w, "npc"),
      height = unit(trait_h, "npc"),
      r = unit(0.010, "npc"),
      gp = gpar(fill = "#111827", col = NA)
    )
    grid.text(
      paste0("Trait\n", trait_id),
      x = unit(trait_x, "npc"),
      y = unit(trait_y + 0.002, "npc"),
      gp = gpar(col = "white", fontsize = 10.0, fontface = "bold", lineheight = 0.92)
    )

    row_top <- 0.775
    row_bottom <- 0.215
    row_step <- if (n_rows <= 1) 0 else (row_top - row_bottom) / (n_rows - 1)
    module_h <- if (n_rows <= 1) 0.135 else min(0.135, max(0.060, row_step * 0.74))

    left_module_x <- 0.195
    left_module_w <- 0.245
    left_program_x <- 0.350
    right_program_x <- 0.650
    right_module_x <- 0.805
    right_module_w <- 0.245

    if (n_left > 0) {
      for (i in seq_len(n_left)) {
        row <- panel_programs[i, , drop = FALSE]
        y <- if (n_rows <= 1) trait_y else row_top - (i - 1) * row_step
        if (i %% 2 == 0) {
          grid.roundrect(
            x = unit(0.250, "npc"),
            y = unit(y, "npc"),
            width = unit(0.375, "npc"),
            height = unit(module_h * 1.12, "npc"),
            r = unit(0.006, "npc"),
            gp = gpar(fill = "#F9FAFB", col = NA)
          )
        }
        hits <- side_hits[side_hits$Program == row$Program & side_hits$side == "program_loading", , drop = FALSE]
        lines <- split_gene_lines(hits)
        draw_gene_module(
          x = left_module_x,
          y = y,
          width = left_module_w,
          height = module_h,
          title = sprintf("%s program genes", row$Program),
          lines = lines,
          gamma_max = gamma_max,
          align = "right",
          fill = "#FFFFFF"
        )
        grid.segments(
          x0 = unit(left_module_x + left_module_w / 2 + 0.006, "npc"),
          x1 = unit(left_program_x - 0.038, "npc"),
          y0 = unit(y, "npc"),
          y1 = unit(y, "npc"),
          gp = gpar(col = "#D1D5DB", lwd = 0.8)
        )
        draw_program_box(
          x = left_program_x,
          y = y,
          label = row$Program,
          score_label = sprintf("B %s", format_signed_score(row$program_score)),
          fill = evidence_fill[as.character(row$color)]
        )
        draw_trait_arrow(
          x0 = left_program_x + 0.038,
          y0 = y,
          x1 = trait_left - 0.004,
          y1 = trait_y,
          direction = row$program_score,
          col = "#4B5563",
          lwd = 1.6 + 1.1 * min(1, abs(as.numeric(row$program_score)) / 4)
        )
      }
    }

    if (n_right > 0) {
      for (i in seq_len(n_right)) {
        row <- panel_regulators[i, , drop = FALSE]
        y <- if (n_rows <= 1) trait_y else row_top - (i - 1) * row_step
        if (i %% 2 == 0) {
          grid.roundrect(
            x = unit(0.750, "npc"),
            y = unit(y, "npc"),
            width = unit(0.375, "npc"),
            height = unit(module_h * 1.12, "npc"),
            r = unit(0.006, "npc"),
            gp = gpar(fill = "#F9FAFB", col = NA)
          )
        }
        hits <- side_hits[side_hits$Program == row$Program & side_hits$side == "regulator", , drop = FALSE]
        lines <- split_gene_lines(hits)
        draw_program_box(
          x = right_program_x,
          y = y,
          label = row$Program,
          score_label = sprintf("R %s", format_signed_score(row$regulator_score)),
          fill = evidence_fill[as.character(row$color)]
        )
        draw_trait_arrow(
          x0 = right_program_x - 0.038,
          y0 = y,
          x1 = trait_right + 0.004,
          y1 = trait_y,
          direction = row$regulator_score,
          col = "#4B5563",
          lwd = 1.6 + 1.1 * min(1, abs(as.numeric(row$regulator_score)) / 4)
        )
        grid.segments(
          x0 = unit(right_program_x + 0.038, "npc"),
          x1 = unit(right_module_x - right_module_w / 2 - 0.006, "npc"),
          y0 = unit(y, "npc"),
          y1 = unit(y, "npc"),
          gp = gpar(col = "#D1D5DB", lwd = 0.8)
        )
        draw_gene_module(
          x = right_module_x,
          y = y,
          width = right_module_w,
          height = module_h,
          title = sprintf("%s regulators", row$Program),
          lines = lines,
          gamma_max = gamma_max,
          align = "left",
          fill = "#FFFFFF"
        )
      }
    }

    legend_y <- 0.095
    grid.points(
      x = unit(c(0.150, 0.255, 0.355), "npc"),
      y = unit(rep(legend_y, 3), "npc"),
      pch = 16,
      size = unit(0.018, "npc"),
      gp = gpar(col = c("#2166AC", "#7A808A", "#B2182B"))
    )
    grid.text("negative gamma", x = unit(0.170, "npc"), y = unit(legend_y, "npc"), just = "left", gp = gpar(col = "#2166AC", fontsize = 7.4, fontface = "bold"))
    grid.text("near 0", x = unit(0.275, "npc"), y = unit(legend_y, "npc"), just = "left", gp = gpar(col = "#7A808A", fontsize = 7.4, fontface = "bold"))
    grid.text("positive gamma", x = unit(0.375, "npc"), y = unit(legend_y, "npc"), just = "left", gp = gpar(col = "#B2182B", fontsize = 7.4, fontface = "bold"))
    grid.text(
      "Each program/regulator module is connected to the center trait; arrow head = positive direction, flat cap = negative direction.",
      x = unit(0.555, "npc"),
      y = unit(legend_y, "npc"),
      just = "left",
      gp = gpar(col = "#374151", fontsize = 7.3)
    )
    grid.text(
      sprintf(
        "Filters: |gamma| >= %.3g; program genes from top %d loading genes; regulators at FDR <= %.3g. B = program burden score; R = regulator-burden correlation score.",
        hit_abs_gamma_threshold,
        loading_top_n,
        regulator_fdr_threshold
      ),
      x = unit(0.020, "npc"),
      y = unit(0.040, "npc"),
      just = "left",
      gp = gpar(col = "#4B5563", fontsize = 7.2)
    )
  }

  grDevices::pdf(paste0(plot_prefix, ".pdf"), width = canvas_width, height = canvas_height, useDingbats = FALSE)
  draw_centered_page()
  grDevices::dev.off()

  grDevices::png(paste0(plot_prefix, ".png"), width = canvas_width, height = canvas_height, units = "in", res = 300)
  draw_centered_page()
  grDevices::dev.off()
}

draw_model_panel(
  program_summary = program_summary,
  side_hits = side_hits,
  trait_id = trait_id,
  plot_prefix = plot_prefix,
  max_genes_per_side = max_genes_per_side,
  hit_abs_gamma_threshold = hit_abs_gamma_threshold,
  loading_top_n = loading_top_n,
  regulator_fdr_threshold = regulator_fdr_threshold
)

quit(save = "no")

signed_value_color <- function(values, max_abs) {
  values <- as.numeric(values)
  values[!is.finite(values)] <- 0
  max_abs <- max(as.numeric(max_abs), 1e-12)

  pos_palette <- grDevices::colorRampPalette(c("#F4B9A9", "#B2182B"))(101)
  neg_palette <- grDevices::colorRampPalette(c("#B9D8EA", "#2166AC"))(101)
  intensity <- pmin(1, abs(values) / max_abs)
  idx <- pmax(1, pmin(101, floor(intensity * 100) + 1))

  colors <- rep("#8A8F98", length(values))
  colors[values > 0] <- pos_palette[idx[values > 0]]
  colors[values < 0] <- neg_palette[idx[values < 0]]
  colors
}

cap_nonfinite_scores <- function(values) {
  values <- as.numeric(values)
  finite_values <- values[is.finite(values)]
  cap <- if (length(finite_values) > 0) max(abs(finite_values), 1, na.rm = TRUE) else 1
  values[is.infinite(values)] <- sign(values[is.infinite(values)]) * cap
  values[is.na(values)] <- 0
  values
}

format_signed_score <- function(values) {
  values <- as.numeric(values)
  out <- ifelse(
    is.na(values),
    "NA",
    ifelse(
      is.finite(values),
      sprintf("%+.1f", values),
      sprintf("%sInf", ifelse(values > 0, "+", "-"))
    )
  )
  out
}

program_plot <- program_summary
program_plot$row_index <- seq_len(panel_rows)
row_spacing <- max(1.75, 0.85 + max_genes_per_side * 0.18)
program_plot$y_row <- (panel_rows - program_plot$row_index + 1) * row_spacing
program_plot$program_x <- 0

program_abs <- abs(program_plot$program_score)
regulator_abs <- abs(program_plot$regulator_score)
program_abs[is.na(program_abs)] <- -Inf
regulator_abs[is.na(regulator_abs)] <- -Inf
use_program_score <- program_abs >= regulator_abs
use_program_score[!is.finite(program_abs) & !is.finite(regulator_abs)] <- TRUE
program_plot$dominant_source <- ifelse(use_program_score, "B", "R")
program_plot$program_effect_raw <- ifelse(use_program_score, program_plot$program_score, program_plot$regulator_score)
program_plot$program_effect_score <- cap_nonfinite_scores(program_plot$program_effect_raw)
program_effect_max <- max(abs(program_plot$program_effect_score), 1, na.rm = TRUE)
program_plot$program_direction <- sign(program_plot$program_effect_score)
program_plot$effect_color <- "#374151"
program_plot$effect_alpha <- 0.32 + 0.58 * pmin(1, abs(program_plot$program_effect_score) / program_effect_max)
program_plot$effect_width <- 0.45 + 1.45 * pmin(1, abs(program_plot$program_effect_score) / program_effect_max)
program_plot$program_node_label <- program_plot$Program
program_plot$score_label <- sprintf(
  "B:%s  R:%s",
  format_signed_score(program_plot$program_score),
  format_signed_score(program_plot$regulator_score)
)
program_plot$effect_y <- program_plot$y_row + row_spacing * 0.20
program_plot$effect_x0 <- -0.38
program_plot$effect_x1 <- 0.38

max_side_count <- max(table(paste(side_hits$Program, side_hits$side, sep = "::")), 1)
gene_y_gap <- min(0.22, (row_spacing * 0.72) / max(1, max_side_count - 1))
gene_plot <- merge(
  side_hits,
  program_plot[, c("Program", "y_row")],
  by = "Program",
  all.x = TRUE,
  sort = FALSE
)
gene_plot <- gene_plot[!is.na(gene_plot$y_row), , drop = FALSE]
gene_plot$side_label <- ifelse(gene_plot$side == "program_loading", "Program loading genes", "Regulator genes")
gene_plot$side_count <- ave(
  seq_len(nrow(gene_plot)),
  gene_plot$Program,
  gene_plot$side,
  FUN = function(x) length(x)
)
gene_plot$display_rank <- as.numeric(gene_plot$display_rank)
gene_plot$y_gene <- gene_plot$y_row + (gene_plot$display_rank - ((gene_plot$side_count + 1) / 2)) * gene_y_gap
gene_plot$gene_point_x <- ifelse(gene_plot$side == "program_loading", -3.45, 3.45)
gene_plot$label_x <- ifelse(gene_plot$side == "program_loading", -3.78, 3.78)
gene_plot$program_edge_x <- ifelse(gene_plot$side == "program_loading", -0.82, 0.82)
gene_plot$label_hjust <- ifelse(gene_plot$side == "program_loading", 1, 0)
gene_plot$post_mean_plot <- ifelse(is.finite(gene_plot$post_mean), gene_plot$post_mean, 0)
gene_gamma_max <- max(abs(gene_plot$post_mean_plot), hit_abs_gamma_threshold, na.rm = TRUE)
gene_plot$gene_color <- signed_value_color(gene_plot$post_mean_plot, gene_gamma_max)
gene_plot$branch_alpha <- 0.18 + 0.24 * pmin(1, abs(gene_plot$post_mean_plot) / gene_gamma_max)
gene_plot$point_size <- 1.65 + 1.15 * pmin(1, abs(gene_plot$post_mean_plot) / gene_gamma_max)

pos_effects <- program_plot[program_plot$program_direction > 0, , drop = FALSE]
neg_effects <- program_plot[program_plot$program_direction < 0, , drop = FALSE]
zero_effects <- program_plot[program_plot$program_direction == 0, , drop = FALSE]
if (nrow(neg_effects) > 0) {
  neg_effects$cap_x <- neg_effects$effect_x1
  neg_effects$cap_y0 <- neg_effects$effect_y - 0.08
  neg_effects$cap_y1 <- neg_effects$effect_y + 0.08
}

y_header <- max(program_plot$y_row) + row_spacing * 0.82
y_trait <- y_header + row_spacing * 0.62
y_legend <- min(program_plot$y_row) - row_spacing * 0.72
plot_y_min <- y_legend - row_spacing * 0.58
plot_y_max <- y_trait + row_spacing * 0.50
plot_x_min <- -5.15
plot_x_max <- 5.15
plot_width <- max(12.8, min(17.0, 11.8 + max_genes_per_side * 0.28))
plot_height <- max(7.8, min(16.5, 3.2 + panel_rows * 0.92))
gene_text_size <- max(2.6, min(3.7, 4.2 - max_genes_per_side * 0.11))
program_text_size <- max(3.0, min(4.0, 4.25 - panel_rows * 0.045))

heading_df <- data.frame(
  x = c(-3.45, 0, 3.45),
  y = y_header,
  label = c(
    sprintf("Program loading genes\n(top %d overlap)", loading_top_n),
    "Selected programs\nB/R signed -log10(P)",
    sprintf("Regulator genes\n(FDR <= %.3g)", regulator_fdr_threshold)
  ),
  stringsAsFactors = FALSE
)
trait_df <- data.frame(x = 0, y = y_trait, label = paste0("Trait: ", trait_id), stringsAsFactors = FALSE)
spine_df <- data.frame(x = 0, y = min(program_plot$y_row) - row_spacing * 0.20, yend = y_trait - row_spacing * 0.32)
row_band_df <- program_plot[program_plot$row_index %% 2 == 0, , drop = FALSE]
row_band_df$xmin <- rep(plot_x_min + 0.16, nrow(row_band_df))
row_band_df$xmax <- rep(plot_x_max - 0.16, nrow(row_band_df))
row_band_df$ymin <- row_band_df$y_row - row_spacing * 0.36
row_band_df$ymax <- row_band_df$y_row + row_spacing * 0.36
gene_color_legend <- data.frame(
  x = c(-4.10, -2.75, -1.60),
  y = y_legend,
  label = c("negative gamma", "near 0", "positive gamma"),
  color = c("#2166AC", "#8A8F98", "#B2182B"),
  stringsAsFactors = FALSE
)
direction_legend <- data.frame(
  x = -0.25,
  y = y_legend,
  label = "program mark: arrow positive, cap negative; (gene) discordant",
  stringsAsFactors = FALSE
)
evidence_label <- data.frame(
  x = 2.55,
  y = y_legend,
  label = "evidence:",
  stringsAsFactors = FALSE
)
evidence_legend <- data.frame(
  x = c(3.12, 3.66, 4.20, 4.74),
  y = y_legend,
  label = c("none", "B", "R", "B+R"),
  color = factor(
    c("other", "program_enriched", "regulator_enriched", "both_enriched"),
    levels = c("other", "program_enriched", "regulator_enriched", "both_enriched")
  ),
  stringsAsFactors = FALSE
)

subtitle_text <- sprintf(
  "%d selected programs; %d genes shown (%d loading, %d regulator). Gene color encodes signed gamma and intensity.",
  panel_rows,
  nrow(gene_plot),
  sum(gene_plot$side == "program_loading"),
  sum(gene_plot$side == "regulator")
)
caption_text <- sprintf(
  "B = program burden score; R = regulator-burden correlation score. Program fill uses Bonferroni evidence: gray = neither, B = program burden, R = regulator burden, B+R = both. Gene color shows signed gamma. Filters: |gamma| >= %.3g, top %d loading genes, regulator FDR <= %.3g.",
  hit_abs_gamma_threshold,
  loading_top_n,
  regulator_fdr_threshold
)

g <- ggplot()
g <- g + geom_rect(
  data = row_band_df,
  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  fill = "#F9FAFB",
  color = NA
)
g <- g + geom_segment(
  data = spine_df,
  aes(x = x, y = y, xend = x, yend = yend),
  linewidth = 0.45,
  color = "#D1D5DB"
)
g <- g + geom_segment(
  data = gene_plot,
  aes(
    x = gene_point_x,
    y = y_gene,
    xend = program_edge_x,
    yend = y_row,
    alpha = branch_alpha,
    linetype = side_label
  ),
  linewidth = 0.20,
  color = "#CBD5E1"
)
g <- g + geom_segment(
  data = zero_effects,
  aes(x = effect_x0, y = effect_y, xend = effect_x1, yend = effect_y, color = effect_color, linewidth = effect_width, alpha = effect_alpha)
)
g <- g + geom_segment(
  data = neg_effects,
  aes(x = effect_x0, y = effect_y, xend = effect_x1, yend = effect_y, color = effect_color, linewidth = effect_width, alpha = effect_alpha)
)
if (nrow(neg_effects) > 0) {
  g <- g + geom_segment(
    data = neg_effects,
    aes(x = cap_x, y = cap_y0, xend = cap_x, yend = cap_y1, color = effect_color, linewidth = effect_width, alpha = effect_alpha)
  )
}
g <- g + geom_segment(
  data = pos_effects,
  aes(x = effect_x0, y = effect_y, xend = effect_x1, yend = effect_y, color = effect_color, linewidth = effect_width, alpha = effect_alpha),
  arrow = grid::arrow(length = grid::unit(0.10, "inches"), type = "closed")
)
g <- g + geom_label(
  data = program_plot,
  aes(x = program_x, y = y_row, label = program_node_label, fill = color),
  label.size = 0.18,
  label.padding = grid::unit(0.18, "lines"),
  lineheight = 0.90,
  size = program_text_size,
  fontface = "bold",
  color = "#111827"
)
g <- g + geom_text(
  data = program_plot,
  aes(x = program_x, y = y_row - row_spacing * 0.25, label = score_label),
  size = max(2.25, program_text_size * 0.66),
  color = "#4B5563",
  fontface = "bold"
)
g <- g + geom_point(
  data = gene_plot,
  aes(x = gene_point_x, y = y_gene, shape = side_label, size = point_size, color = gene_color),
  stroke = 0,
  alpha = 0.78
)
g <- g + geom_text(
  data = gene_plot,
  aes(
    x = label_x,
    y = y_gene,
    label = gene_label,
    hjust = label_hjust,
    color = gene_color,
    fontface = ifelse(is_discordant, "bold", "plain")
  ),
  size = gene_text_size,
  lineheight = 0.9
)
g <- g + geom_label(
  data = trait_df,
  aes(x = x, y = y, label = label),
  fill = "#111827",
  color = "white",
  label.size = 0,
  label.padding = grid::unit(0.25, "lines"),
  size = 4.2,
  fontface = "bold"
)
g <- g + geom_text(
  data = heading_df,
  aes(x = x, y = y, label = label),
  size = 3.35,
  lineheight = 0.92,
  fontface = "bold",
  color = "#111827"
)
g <- g + geom_point(
  data = gene_color_legend,
  aes(x = x, y = y, color = color),
  size = 3.1
)
g <- g + geom_text(
  data = gene_color_legend,
  aes(x = x + 0.12, y = y, label = label, color = color),
  hjust = 0,
  size = 2.75,
  fontface = "bold"
)
g <- g + geom_label(
  data = direction_legend,
  aes(x = x, y = y, label = label),
  hjust = 0,
  vjust = 0.5,
  fill = "white",
  color = "#374151",
  label.size = 0,
  label.padding = grid::unit(0.12, "lines"),
  lineheight = 0.95,
  size = 2.45
)
g <- g + geom_text(
  data = evidence_label,
  aes(x = x, y = y, label = label),
  hjust = 0,
  vjust = 0.5,
  size = 2.45,
  fontface = "bold",
  color = "#374151"
)
g <- g + geom_label(
  data = evidence_legend,
  aes(x = x, y = y, label = label, fill = color),
  label.size = 0.12,
  label.padding = grid::unit(0.12, "lines"),
  size = 2.25,
  fontface = "bold",
  color = "#111827"
)
g <- g + scale_color_identity()
g <- g + scale_alpha_identity()
g <- g + scale_linewidth_identity()
g <- g + scale_fill_manual(
  values = c(
    "other" = "#F3F4F6",
    "program_enriched" = "#F9D48B",
    "regulator_enriched" = "#A9C7E6",
    "both_enriched" = "#A9D6A4"
  ),
  drop = FALSE,
  name = "Program evidence"
)
g <- g + scale_shape_manual(
  values = c("Program loading genes" = 16, "Regulator genes" = 17),
  name = "Gene source"
)
g <- g + scale_linetype_manual(
  values = c("Program loading genes" = "solid", "Regulator genes" = "22"),
  name = "Gene source"
)
g <- g + scale_size_identity()
g <- g + coord_cartesian(xlim = c(plot_x_min, plot_x_max), ylim = c(plot_y_min, plot_y_max), clip = "off")
g <- g + labs(
  title = "Programs selected for modeling trait associations",
  subtitle = subtitle_text,
  caption = caption_text
)
g <- g + theme_void(base_family = "Helvetica")
g <- g + theme(
  plot.background = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA),
  legend.position = "none",
  legend.box = "vertical",
  legend.title = element_text(size = 9, face = "bold", color = "#111827"),
  legend.text = element_text(size = 8, color = "#374151"),
  plot.title = element_text(size = 17, face = "bold", color = "#111827", hjust = 0.5),
  plot.subtitle = element_text(size = 9.5, color = "#374151", hjust = 0.5, margin = margin(t = 4, b = 10)),
  plot.caption = element_text(size = 8.2, color = "#4B5563", hjust = 0, margin = margin(t = 10)),
  plot.margin = margin(18, 30, 18, 30)
)

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = plot_width, height = plot_height)
ggsave(plot = g, filename = paste0(plot_prefix, ".png"), width = plot_width, height = plot_height, dpi = 300)
