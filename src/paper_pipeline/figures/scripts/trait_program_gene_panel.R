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

head_by_program <- function(df, max_rows) {
  if (is.null(df) || nrow(df) == 0) {
    return(df[0, , drop = FALSE])
  }
  groups <- split(df, df$Program, drop = TRUE)
  if (length(groups) == 0) {
    return(df[0, , drop = FALSE])
  }
  out <- do.call(rbind, lapply(groups, function(group_df) head(group_df, max_rows)))
  row.names(out) <- NULL
  out
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

trait_file <- basename(posterior_path)
program_summary <- read_program_regulator_summary(program_association_dir, trait_file, k, character())
program_summary$Program <- as.character(program_summary$Program)
program_summary$program_sig <- program_summary$MEANgamma_top100_shet_adjusted_P < 0.05 / nrow(program_summary)
program_summary$regulator_sig <- program_summary$P_withShet < 0.05 / nrow(program_summary)
program_summary$priority_tier <- ifelse(
  program_summary$program_sig & program_summary$regulator_sig, 1,
  ifelse(program_summary$program_sig | program_summary$regulator_sig, 2, 3)
)
program_summary$priority_score <- pmax(abs(program_summary$program_score), abs(program_summary$regulator_score))
program_summary <- program_summary[order(program_summary$priority_tier, -program_summary$priority_score, program_summary$Program), , drop = FALSE]
program_summary <- program_summary[
  program_summary$priority_tier < 3 | program_summary$priority_score >= min_abs_score,
  ,
  drop = FALSE
]
if (nrow(program_summary) == 0) {
  program_summary <- read_program_regulator_summary(program_association_dir, trait_file, k, character())
  program_summary$Program <- as.character(program_summary$Program)
  program_summary$priority_score <- pmax(abs(program_summary$program_score), abs(program_summary$regulator_score))
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
  loading_hits <- head_by_program(loading_hits, max_genes_per_side)
} else {
  loading_hits <- loading_hits[0, , drop = FALSE]
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
  regulator_hits <- head_by_program(regulator_hits, max_genes_per_side)
} else {
  regulator_hits <- regulator_hits[0, , drop = FALSE]
}

side_hit_cols <- c("Program", "side", "gene", "ensg", "post_mean", "abs_gamma", "gamma_sign", "membership_score", "rank_within_side")
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

gene_point_size <- max(2.8, min(5.0, 6.0 - max_genes_per_side * 0.18))
gene_text_size <- max(3.2, min(5.2, 5.8 - max_genes_per_side * 0.22))
plot_height <- max(7, panel_rows * (1.0 + max_genes_per_side * 0.25))

g <- ggplot()
g <- g + geom_segment(
  data = program_summary,
  aes(x = -0.78, xend = 0.78, y = y_center, yend = y_center),
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
g <- g + geom_point(
  data = side_hits,
  aes(x = x, y = y_global, color = gamma_sign),
  size = gene_point_size
)
g <- g + geom_text(
  data = side_hits,
  aes(
    x = ifelse(side_hits$side == "program_loading", -0.94, 0.94),
    y = y_global,
    label = gene_label,
    hjust = ifelse(side_hits$side == "program_loading", 1, 0),
    color = gamma_sign
  ),
  size = gene_text_size
)
g <- g + annotate("text", x = -1, y = max(side_hits$y_global) + 2, label = "Program loading genes", hjust = 0.5, fontface = "bold", size = gene_text_size * 1.05)
g <- g + annotate("text", x = 1, y = max(side_hits$y_global) + 2, label = "Program regulators", hjust = 0.5, fontface = "bold", size = gene_text_size * 1.05)
g <- g + annotate("text", x = 0, y = max(side_hits$y_global) + 2, label = trait_id, hjust = 0.5, fontface = "bold", size = gene_text_size * 1.1)
g <- g + scale_color_manual(values = c("positive" = "#B40426", "negative" = "#3B4CC0", "zero" = "grey55"))
g <- g + scale_fill_manual(
  values = c(
    "other" = "#ECECEC",
    "program_enriched" = "#F9D48B",
    "regulator_enriched" = "#A9C7E6",
    "both_enriched" = "#A9D6A4"
  ),
  drop = FALSE
)
g <- g + coord_cartesian(xlim = c(-1.28, 1.28), ylim = c(0, max(side_hits$y_global) + 3), clip = "off")
g <- g + theme_void(base_family = "Helvetica")
g <- g + theme(
  legend.position = "none",
  plot.margin = margin(18, 28, 18, 28)
)

ggsave(plot = g, filename = paste0(plot_prefix, ".pdf"), width = 13, height = plot_height)
ggsave(plot = g, filename = paste0(plot_prefix, ".png"), width = 13, height = plot_height, dpi = 300)
