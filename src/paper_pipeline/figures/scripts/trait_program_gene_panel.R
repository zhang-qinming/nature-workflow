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
max_programs <- if (length(args) >= 10) as.numeric(args[10]) else 8
max_genes_per_side <- if (length(args) >= 11) as.numeric(args[11]) else 8
hit_abs_gamma_threshold <- if (length(args) >= 12) as.numeric(args[12]) else 0.1
loading_top_n <- if (length(args) >= 13) as.numeric(args[13]) else 200
regulator_fdr_threshold <- if (length(args) >= 14) as.numeric(args[14]) else 0.05
min_abs_score <- if (length(args) >= 15) as.numeric(args[15]) else 1.3
render_plot <- if (length(args) >= 16) as.character(args[16]) else "1"
association_trait_file <- if (length(args) >= 17) as.character(args[17]) else ""

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
    predicted_sign = numeric(),
    post_mean_sign = numeric(),
    is_discordant = logical(),
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
program_summary <- read_program_regulator_summary(program_association_dir, trait_file, k, character())
program_summary <- annotate_program_summary(program_summary)
program_summary <- program_summary[order(program_summary$priority_tier, -program_summary$priority_score, program_summary$Program), , drop = FALSE]
program_summary <- program_summary[
  program_summary$priority_tier < 3 | program_summary$priority_score >= min_abs_score,
  ,
  drop = FALSE
]
if (nrow(program_summary) == 0) {
  program_summary <- read_program_regulator_summary(program_association_dir, trait_file, k, character())
  program_summary <- annotate_program_summary(program_summary)
  program_summary <- program_summary[order(-program_summary$priority_score, program_summary$Program), , drop = FALSE]
}
program_summary <- head(program_summary, max_programs)
selected_programs <- as.character(program_summary$Program)

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

trait_hits <- posterior_df[!is.na(posterior_df$abs_gamma) & posterior_df$abs_gamma >= hit_abs_gamma_threshold, , drop = FALSE]
if (nrow(trait_hits) == 0) {
  trait_hits <- posterior_df[order(-posterior_df$abs_gamma), , drop = FALSE]
  trait_hits <- head(trait_hits, min(30, nrow(trait_hits)))
}

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

if (has_overlap) {
  program_mean_sign <- sign(program_summary$MEANgamma_top100 - program_summary$shet_adjusted_random_mean)
  names(program_mean_sign) <- program_summary$Program

  regulator_sign_lookup <- if (nrow(regulator_df) > 0) {
    setNames(sign(regulator_df$beta), paste(regulator_df$Program, regulator_df$gene, sep = "::"))
  } else {
    character()
  }

  side_hits$predicted_sign <- 0
  is_loading <- side_hits$side == "program_loading"
  side_hits$predicted_sign[is_loading] <- unname(program_mean_sign[side_hits$Program[is_loading]])
  if (any(!is_loading)) {
    reg_keys <- paste(side_hits$Program[!is_loading], side_hits$gene[!is_loading], sep = "::")
    side_hits$predicted_sign[!is_loading] <- unname(regulator_sign_lookup[reg_keys])
  }
  side_hits$predicted_sign[is.na(side_hits$predicted_sign)] <- 0
  side_hits$post_mean_sign <- sign(side_hits$post_mean)
  side_hits$is_discordant <- side_hits$predicted_sign != 0 & side_hits$post_mean_sign != 0 & side_hits$predicted_sign == -side_hits$post_mean_sign

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
