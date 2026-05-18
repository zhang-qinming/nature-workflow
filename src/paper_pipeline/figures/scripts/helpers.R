trait_stem <- function(trait_file) {
  stem <- basename(as.character(trait_file))
  sub("\\.per_gene_estimates\\.tsv$", "", stem)
}

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

parse_program_vector <- function(raw) {
  if (is.null(raw) || !nzchar(raw)) {
    return(character())
  }
  paste0("P", trimws(strsplit(as.character(raw), ",", fixed = TRUE)[[1]]))
}

friendly_geneset_name <- function(name) {
  mapping <- c(
    HALLMARK_HEME_METABOLISM = "Heme metabolism",
    Hematopoiesisgenes = "Hematopoiesis",
    mitotic_cell_cycle = "Mitotic cell cycle",
    positive_macromolecule_synthesis = "Macromolecule biosynthesis, positive reg"
  )
  ifelse(name %in% names(mapping), unname(mapping[name]), name)
}

resolve_geneset_path <- function(geneset_dir, geneset_name) {
  direct <- file.path(geneset_dir, geneset_name)
  with_txt <- file.path(geneset_dir, paste0(geneset_name, ".txt"))
  if (file.exists(direct)) {
    return(direct)
  }
  if (file.exists(with_txt)) {
    return(with_txt)
  }
  stop(sprintf("Geneset file not found for '%s' in %s", geneset_name, geneset_dir))
}

read_geneset_members <- function(geneset_dir, geneset_name) {
  geneset_path <- resolve_geneset_path(geneset_dir, geneset_name)
  unique(read.table(geneset_path, header = FALSE, stringsAsFactors = FALSE)[, 1])
}

read_gene_map <- function(gene_map_path) {
  gene_map <- utils::read.table(gene_map_path, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(gene_map) < 2) {
    stop(sprintf("Gene map at %s must have at least 2 columns", gene_map_path))
  }
  gene_map <- gene_map[, 1:2]
  colnames(gene_map) <- c("ensg", "gene")
  gene_map
}

label_from_ensg <- function(ensg, gene_map_path) {
  gene_map <- read_gene_map(gene_map_path)
  lookup <- stats::setNames(gene_map$gene, gene_map$ensg)
  labels <- lookup[ensg]
  labels[is.na(labels)] <- ensg[is.na(labels)]
  unname(labels)
}

parse_variant_table <- function(gwas) {
  # auto-detect column names (case-insensitive match)
  nms <- names(gwas)
  find_col <- function(candidates) {
    for (col in candidates) {
      idx <- which(tolower(nms) == tolower(col))
      if (length(idx) == 1) return(nms[idx])
    }
    return(NULL)
  }

  # detect variant identifier
  variant_col <- find_col(c("variant", "SNP", "snp", "rsid", "MarkerName"))
  chrom_col   <- find_col(c("CHROM", "chr", "CHR", "chromosome"))
  pos_col     <- find_col(c("POS", "pos", "BP", "bp", "position"))
  pval_col    <- find_col(c("pval", "p_value", "p.value", "P", "p", "pvalue"))

  if (is.null(pval_col)) {
    stop(sprintf("GWAS file missing p-value column. Available: %s",
                 paste(nms, collapse = ", ")))
  }

  gwas_f <- gwas
  if ("low_confidence_variant" %in% nms) {
    gwas_f <- gwas_f[gwas_f$low_confidence_variant == FALSE, , drop = FALSE]
  }

  # parse chromosome and position
  if (!is.null(variant_col) && grepl(":", as.character(gwas_f[[variant_col]][1]))) {
    # format: CHR:POS:REF:ALT
    variant_parts <- strsplit(as.character(gwas_f[[variant_col]]), ":", fixed = TRUE)
    chrom <- sub("^chr", "", vapply(variant_parts, function(x) x[[1]], character(1)))
    pos   <- suppressWarnings(as.numeric(vapply(variant_parts, function(x) x[[2]], character(1))))
  } else if (!is.null(chrom_col) && !is.null(pos_col)) {
    # separate CHR + POS columns
    chrom <- sub("^chr", "", as.character(gwas_f[[chrom_col]]))
    pos   <- suppressWarnings(as.numeric(gwas_f[[pos_col]]))
  } else if (!is.null(variant_col)) {
    # rsID or other non-parseable variant — try chr/pos from column names
    stop(sprintf("GWAS file has variant/SNP column but cannot parse CHR:POS. Need CHR+POS columns. Available: %s",
                 paste(nms, collapse = ", ")))
  } else {
    stop(sprintf("GWAS file format not recognized. Need: variant (CHR:POS:REF:ALT) or CHR+POS columns. Available: %s",
                 paste(nms, collapse = ", ")))
  }

  # strip leading zeros from chromosome (e.g. "01" → "1")
  chrom <- sub("^0+", "", chrom)

  df <- data.frame(
    CHROM = chrom,
    POS = pos,
    P = as.numeric(gwas_f[[pval_col]]),
    stringsAsFactors = FALSE
  )
  n_before <- nrow(df)
  df <- df[!is.na(df$POS) & !is.na(df$P), , drop = FALSE]
  n_after_na <- nrow(df)
  # keep autosomes only (1-22)
  df <- df[df$CHROM %in% as.character(seq_len(22)), , drop = FALSE]
  n_final <- nrow(df)

  message(sprintf("GWAS parsing: %d variants loaded, %d after NA filter, %d autosomes kept (dropped %d non-autosomal)",
                  n_before, n_after_na, n_final, n_after_na - n_final))

  if (n_final == 0) {
    if (n_after_na == 0) {
      stop(sprintf("No valid variants after NA filter (%d raw). Check POS and P columns for missing data.", n_before))
    } else {
      stop(sprintf("No autosomal variants found (%d non-autosomal). Check chromosome format. Unique chroms: %s",
                   n_after_na, paste(unique(chrom[!is.na(chrom)]), collapse = ", ")))
    }
  }
  df
}

read_gene_annotation_table <- function(gene_annotation_path) {
  suppressWarnings({
    gene_annotation <- tryCatch(
      data.table::fread(gene_annotation_path, data.table = FALSE),
      error = function(...) NULL
    )
  })

  required_columns <- c("gene", "chr", "start", "end")
  if (!is.null(gene_annotation) && all(required_columns %in% names(gene_annotation))) {
    return(gene_annotation[, required_columns])
  }

  raw_lines <- readLines(gene_annotation_path, warn = FALSE)
  raw_lines <- raw_lines[!grepl("^#", raw_lines)]
  if (length(raw_lines) == 0) {
    stop(sprintf("Gene annotation %s is empty", gene_annotation_path))
  }

  con <- textConnection(raw_lines)
  on.exit(close(con), add = TRUE)
  gtf <- utils::read.delim(
    con,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = ""
  )
  if (ncol(gtf) < 9) {
    stop(sprintf("Gene annotation %s must be either a table with gene/chr/start/end or a valid GTF", gene_annotation_path))
  }
  colnames(gtf)[1:9] <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  gtf <- gtf[gtf$feature == "gene", , drop = FALSE]
  if (nrow(gtf) == 0) {
    stop(sprintf("No gene features were found in %s", gene_annotation_path))
  }

  gene_name <- rep(NA_character_, nrow(gtf))
  gene_matches <- regmatches(gtf$attribute, regexpr('gene_name "([^"]+)"', gtf$attribute, perl = TRUE))
  gene_name[gene_matches != ""] <- sub('gene_name "([^"]+)"', "\\1", gene_matches[gene_matches != ""])
  gene_id_matches <- regmatches(gtf$attribute, regexpr('gene_id "([^"]+)"', gtf$attribute, perl = TRUE))
  gene_name[is.na(gene_name) & gene_id_matches != ""] <- sub('gene_id "([^"]+)"', "\\1", gene_id_matches[is.na(gene_name) & gene_id_matches != ""])

  annotation_df <- data.frame(
    gene = gene_name,
    chr = gtf$chr,
    start = as.numeric(gtf$start),
    end = as.numeric(gtf$end),
    stringsAsFactors = FALSE
  )
  annotation_df <- annotation_df[!is.na(annotation_df$gene) & !duplicated(annotation_df$gene), , drop = FALSE]
  annotation_df
}

add_cumulative_coordinates <- function(gwas_df) {
  gwas_df <- gwas_df[order(as.numeric(gwas_df$CHROM), gwas_df$POS), , drop = FALSE]
  gwas_df$COORD <- 0

  offset <- 0
  centers <- numeric()
  labels <- character()
  for (chr in seq_len(22)) {
    idx <- gwas_df$CHROM == as.character(chr)
    if (!any(idx)) {
      next
    }
    gwas_df$COORD[idx] <- gwas_df$POS[idx] + offset
    chr_coords <- gwas_df$COORD[idx]
    centers <- c(centers, mean(range(chr_coords)))
    labels <- c(labels, as.character(chr))
    offset <- max(chr_coords)
  }

  list(
    data = gwas_df,
    axis = data.frame(CHROM = labels, center = centers, stringsAsFactors = FALSE)
  )
}

read_program_regulator_summary <- function(association_dir, trait_file, k, label_programs) {
  programs_path <- file.path(association_dir, paste0("programs_enrichment_K", k, "_", trait_file))
  regulators_path <- file.path(association_dir, paste0("regulators_enrichment_K", k, "_", trait_file))

  programs <- data.table::fread(programs_path, data.table = FALSE)
  regulators <- data.table::fread(regulators_path, data.table = FALSE)
  df <- merge(programs, regulators, by = "Program", all = FALSE)

  bonferroni_cutoff <- 0.05 / nrow(df)
  df$program_score <- sign(df$MEANgamma_top100 - df$shet_adjusted_random_mean) * (-log10(df$MEANgamma_top100_shet_adjusted_P))
  df$regulator_score <- sign(df$beta_withShet) * (-log10(df$P_withShet))
  df$label <- df$Program
  df$label[
    df$MEANgamma_top100_shet_adjusted_P > bonferroni_cutoff &
      df$P_withShet > bonferroni_cutoff &
      !is.element(df$Program, label_programs)
  ] <- ""

  df$color <- ifelse(
    df$MEANgamma_top100_shet_adjusted_P < bonferroni_cutoff & df$P_withShet < bonferroni_cutoff,
    "both_enriched",
    ifelse(
      df$MEANgamma_top100_shet_adjusted_P < bonferroni_cutoff,
      "program_enriched",
      ifelse(df$P_withShet < bonferroni_cutoff, "regulator_enriched", "other")
    )
  )
  df$color <- factor(df$color, levels = c("other", "program_enriched", "regulator_enriched", "both_enriched"))
  df
}

read_regulation_matrices <- function(regulation_dir, k) {
  beta_df <- data.frame()
  p_df <- data.frame()

  for (i in seq_len(k)) {
    program_path <- file.path(regulation_dir, paste0("K", k, "_program", i, "_perturb_effects.txt"))
    tmp <- data.table::fread(program_path, data.table = FALSE)

    tmp_beta <- tmp[, 1:2]
    tmp_p <- tmp[, 1:3]

    colnames(tmp_beta) <- c("GENE", paste0("P", i))
    colnames(tmp_p) <- c("GENE", paste0("beta_", i), paste0("P", i))
    tmp_p <- tmp_p[, c("GENE", paste0("P", i))]

    if (i == 1) {
      beta_df <- tmp_beta
      p_df <- tmp_p
    } else {
      beta_df <- merge(beta_df, tmp_beta, by = "GENE", all = FALSE)
      p_df <- merge(p_df, tmp_p, by = "GENE", all = FALSE)
    }
  }

  list(beta = beta_df, p = p_df)
}
