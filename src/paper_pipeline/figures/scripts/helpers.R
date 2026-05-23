trait_stem <- function(trait_file) {
  stem <- basename(as.character(trait_file))
  sub("\\.per_gene_estimates\\.tsv$", "", stem)
}

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

normalize_program_ids <- function(programs) {
  raw <- trimws(as.character(programs))
  raw[is.na(raw)] <- ""
  stripped <- sub("^P", "", raw, ignore.case = TRUE)
  numeric_program <- grepl("^[0-9]+$", stripped)
  raw[numeric_program] <- paste0("P", as.integer(stripped[numeric_program]))
  raw
}

parse_program_vector <- function(raw) {
  if (is.null(raw) || !nzchar(raw)) {
    return(character())
  }
  tokens <- trimws(strsplit(as.character(raw), ",", fixed = TRUE)[[1]])
  normalize_program_ids(tokens[nzchar(tokens)])
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

normalize_ensg_ids <- function(ids) {
  ids <- trimws(as.character(ids))
  ids[is.na(ids)] <- ""
  sub("\\..*$", "", ids)
}

read_gene_map <- function(gene_map_path) {
  gene_map <- utils::read.table(gene_map_path, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(gene_map) < 2) {
    stop(sprintf("Gene map at %s must have at least 2 columns", gene_map_path))
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

build_gene_id_lookups <- function(gene_map_path) {
  gene_map <- read_gene_map(gene_map_path)
  gene_lookup <- stats::setNames(gene_map$gene, gene_map$ensg)
  ensg_lookup <- stats::setNames(gene_map$ensg, gene_map$gene)

  resolve_to_ensg <- function(ids) {
    raw_ids <- trimws(as.character(ids))
    normalized_ids <- normalize_ensg_ids(raw_ids)
    resolved <- ifelse(grepl("^ENSG", normalized_ids), normalized_ids, unname(ensg_lookup[raw_ids]))
    resolved <- normalize_ensg_ids(resolved)
    resolved[!nzchar(resolved)] <- NA_character_
    unname(resolved)
  }

  list(
    gene_map = gene_map,
    gene_lookup = gene_lookup,
    ensg_lookup = ensg_lookup,
    resolve_to_ensg = resolve_to_ensg
  )
}

read_program_membership_lookup <- function(spectra_path, gene_map_path, k = 60, top_n = 100) {
  spectra_raw <- data.table::fread(spectra_path, data.table = FALSE)
  if (ncol(spectra_raw) < 1) {
    stop(sprintf("cNMF spectra at %s is empty", spectra_path))
  }

  if (is.character(spectra_raw[[1]])) {
    row.names(spectra_raw) <- spectra_raw[[1]]
    spectra_raw <- spectra_raw[, -1, drop = FALSE]
  }
  spectra_mat <- t(as.matrix(spectra_raw))
  if (ncol(spectra_mat) == 0) {
    return(character())
  }

  program_count <- min(as.integer(k), ncol(spectra_mat))
  colnames(spectra_mat) <- paste0("P", seq_len(ncol(spectra_mat)))

  gene_map <- read_gene_map(gene_map_path)
  gene_lookup <- stats::setNames(gene_map$gene, gene_map$ensg)
  gene_symbols <- unname(gene_lookup[row.names(spectra_mat)])
  missing_gene <- is.na(gene_symbols) | !nzchar(gene_symbols)
  gene_symbols[missing_gene] <- row.names(spectra_mat)[missing_gene]

  membership <- list()
  for (program_index in seq_len(program_count)) {
    program_name <- paste0("P", program_index)
    program_df <- data.frame(
      gene = gene_symbols,
      weight = as.numeric(spectra_mat[, program_index]),
      stringsAsFactors = FALSE
    )
    program_df <- program_df[order(program_df$weight, decreasing = TRUE, na.last = NA), , drop = FALSE]
    program_df <- program_df[nzchar(program_df$gene) & !duplicated(program_df$gene), , drop = FALSE]
    if (nrow(program_df) == 0) {
      next
    }

    top_genes <- head(program_df$gene, as.integer(top_n))
    for (gene in top_genes) {
      membership[[gene]] <- c(membership[[gene]], program_name)
    }
  }

  if (length(membership) == 0) {
    return(character())
  }

  stats::setNames(
    vapply(membership, function(values) paste(unique(values), collapse = ";"), character(1)),
    names(membership)
  )
}

label_from_ensg <- function(ensg, gene_map_path) {
  gene_map <- read_gene_map(gene_map_path)
  lookup <- stats::setNames(gene_map$gene, gene_map$ensg)
  labels <- lookup[ensg]
  labels[is.na(labels)] <- ensg[is.na(labels)]
  unname(labels)
}

read_program_gene_membership <- function(spectra_path, gene_map_path, k = 60, top_n = 100) {
  spectra_raw <- data.table::fread(spectra_path, data.table = FALSE)
  if (ncol(spectra_raw) < 1) {
    stop(sprintf("cNMF spectra at %s is empty", spectra_path))
  }

  if (is.character(spectra_raw[[1]])) {
    row.names(spectra_raw) <- spectra_raw[[1]]
    spectra_raw <- spectra_raw[, -1, drop = FALSE]
  }
  spectra_mat <- t(as.matrix(spectra_raw))
  if (ncol(spectra_mat) == 0) {
    return(data.frame(program = character(), gene = character(), stringsAsFactors = FALSE))
  }

  program_count <- min(as.integer(k), ncol(spectra_mat))
  colnames(spectra_mat) <- paste0("P", seq_len(ncol(spectra_mat)))

  gene_map <- read_gene_map(gene_map_path)
  gene_lookup <- stats::setNames(gene_map$gene, gene_map$ensg)
  gene_symbols <- unname(gene_lookup[row.names(spectra_mat)])
  missing_gene <- is.na(gene_symbols) | !nzchar(gene_symbols)
  gene_symbols[missing_gene] <- row.names(spectra_mat)[missing_gene]

  rows <- list()
  for (program_index in seq_len(program_count)) {
    program_name <- paste0("P", program_index)
    program_df <- data.frame(
      gene = gene_symbols,
      weight = as.numeric(spectra_mat[, program_index]),
      stringsAsFactors = FALSE
    )
    program_df <- program_df[order(program_df$weight, decreasing = TRUE, na.last = NA), , drop = FALSE]
    program_df <- program_df[nzchar(program_df$gene) & !duplicated(program_df$gene), , drop = FALSE]
    if (nrow(program_df) == 0) {
      next
    }

    top_genes <- head(program_df$gene, as.integer(top_n))
    rows[[length(rows) + 1]] <- data.frame(
      program = program_name,
      gene = top_genes,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) {
    return(data.frame(program = character(), gene = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

annotate_gene_sets <- function(genes, geneset_dir, geneset_names) {
  genes <- as.character(genes)
  output <- rep("", length(genes))
  if (length(geneset_names) == 0 || length(genes) == 0) {
    return(output)
  }

  for (geneset_name in rev(geneset_names)) {
    members <- read_geneset_members(geneset_dir, geneset_name)
    display_name <- friendly_geneset_name(geneset_name)
    mask <- genes %in% members
    if (!any(mask)) {
      next
    }
    output[mask] <- ifelse(
      nchar(output[mask]) == 0,
      display_name,
      paste(output[mask], display_name, sep = ";")
    )
  }

  output
}

annotate_programs <- function(genes, spectra_path, gene_map_path, k = 60, top_n = 100) {
  genes <- as.character(genes)
  output <- rep("", length(genes))
  if (length(genes) == 0) {
    return(output)
  }
  if (!nzchar(spectra_path) || !nzchar(gene_map_path) || !file.exists(spectra_path) || !file.exists(gene_map_path)) {
    return(output)
  }

  lookup <- read_program_membership_lookup(
    spectra_path = spectra_path,
    gene_map_path = gene_map_path,
    k = as.integer(k),
    top_n = as.integer(top_n)
  )
  values <- unname(lookup[genes])
  values[is.na(values)] <- ""
  values
}

annotate_gene_features <- function(
  genes,
  geneset_dir,
  geneset_names = character(),
  spectra_path = "",
  gene_map_path = "",
  k = 60,
  top_n_program_genes = 100
) {
  gene_labels <- as.character(genes)
  genesets <- annotate_gene_sets(gene_labels, geneset_dir, geneset_names)
  programs <- annotate_programs(
    gene_labels,
    spectra_path = spectra_path,
    gene_map_path = gene_map_path,
    k = as.integer(k),
    top_n = as.integer(top_n_program_genes)
  )

  data.frame(
    gene = gene_labels,
    geneset = genesets,
    program = programs,
    stringsAsFactors = FALSE
  )
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
  snp <- if (!is.null(variant_col)) as.character(gwas_f[[variant_col]]) else rep(NA_character_, nrow(gwas_f))

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
  fallback_snp <- paste0(chrom, ":", format(pos, scientific = FALSE, trim = TRUE))
  missing_snp <- is.na(snp) | !nzchar(snp)
  snp[missing_snp] <- fallback_snp[missing_snp]

  df <- data.frame(
    CHROM = chrom,
    POS = pos,
    SNP = snp,
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

  missing_paths <- c(programs_path, regulators_path)[!file.exists(c(programs_path, regulators_path))]
  if (length(missing_paths) > 0) {
    trait_hint <- sub("[._].*$", "", basename(trait_file))
    candidates <- if (dir.exists(association_dir)) {
      list.files(association_dir, pattern = paste0("_K", k, "_"))
    } else {
      character()
    }
    hinted <- candidates[grepl(trait_hint, candidates, fixed = TRUE)]
    shown <- if (length(hinted) > 0) hinted else candidates
    candidate_text <- if (length(shown) > 0) paste(head(shown, 12), collapse = ", ") else "<none>"
    stop(sprintf(
      paste(
        "Missing ProgramLevel association file(s) for trait_file '%s'.",
        "Expected programs file: %s",
        "Expected regulators file: %s",
        "Candidate files in association_dir: %s",
        sep = "\n"
      ),
      trait_file,
      programs_path,
      regulators_path,
      candidate_text
    ))
  }

  programs <- data.table::fread(programs_path, data.table = FALSE)
  regulators <- data.table::fread(regulators_path, data.table = FALSE)
  if (!"Program" %in% names(programs)) {
    stop(sprintf("Program summary %s is missing column: Program", programs_path))
  }
  if (!"Program" %in% names(regulators)) {
    stop(sprintf("Regulator summary %s is missing column: Program", regulators_path))
  }
  programs$Program <- normalize_program_ids(programs$Program)
  regulators$Program <- normalize_program_ids(regulators$Program)
  df <- merge(programs, regulators, by = "Program", all = FALSE)

  if (nrow(df) == 0) {
    stop(sprintf(
      "Program and regulator summaries for trait_file '%s' have no overlapping Program IDs after normalization",
      trait_file
    ))
  }

  required_columns <- c(
    "MEANgamma_top100",
    "shet_adjusted_random_mean",
    "MEANgamma_top100_shet_adjusted_P",
    "beta_withShet",
    "P_withShet"
  )
  missing_columns <- setdiff(required_columns, names(df))
  if (length(missing_columns) > 0) {
    stop(sprintf(
      "ProgramLevel summary for trait_file '%s' is missing columns: %s",
      trait_file,
      paste(missing_columns, collapse = ", ")
    ))
  }

  bonferroni_cutoff <- 0.05 / nrow(df)
  df$program_score <- sign(df$MEANgamma_top100 - df$shet_adjusted_random_mean) * (-log10(df$MEANgamma_top100_shet_adjusted_P))
  df$regulator_score <- sign(df$beta_withShet) * (-log10(df$P_withShet))

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

read_spectra_top_genes <- function(spectra_path, gene_map_path, k, top_n) {
  spectra_raw <- data.table::fread(spectra_path, data.table = FALSE)
  if (ncol(spectra_raw) < 1) {
    stop(sprintf("cNMF spectra at %s is empty", spectra_path))
  }

  if (is.character(spectra_raw[[1]])) {
    row.names(spectra_raw) <- spectra_raw[[1]]
    spectra_raw <- spectra_raw[, -1, drop = FALSE]
  }
  spectra_mat <- t(as.matrix(spectra_raw))
  actual_k <- ncol(spectra_mat)
  colnames(spectra_mat) <- paste0("P", seq_len(actual_k))

  lookups <- build_gene_id_lookups(gene_map_path)
  gene_symbols <- unname(lookups$gene_lookup[normalize_ensg_ids(row.names(spectra_mat))])
  missing_gene <- is.na(gene_symbols) | !nzchar(gene_symbols)
  gene_symbols[missing_gene] <- row.names(spectra_mat)[missing_gene]

  rows <- list()
  for (program_index in seq_len(min(k, actual_k))) {
    program_name <- paste0("P", program_index)
    tmp <- data.frame(
      gene = gene_symbols,
      weight = as.numeric(spectra_mat[, program_index]),
      stringsAsFactors = FALSE
    )
    tmp <- tmp[order(tmp$weight, decreasing = TRUE, na.last = NA), , drop = FALSE]
    tmp <- tmp[nzchar(tmp$gene) & !duplicated(tmp$gene), , drop = FALSE]
    if (nrow(tmp) == 0) {
      next
    }
    tmp <- head(tmp, top_n)
    tmp$Program <- program_name
    tmp$rank_within_side <- seq_len(nrow(tmp))
    rows[[length(rows) + 1]] <- tmp
  }

  if (length(rows) == 0) {
    return(data.frame(
      Program = character(),
      gene = character(),
      weight = numeric(),
      rank_within_side = integer(),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, rows)
}

read_regulator_long <- function(regulation_dir, k) {
  rows <- list()
  for (program_index in seq_len(k)) {
    program_path <- file.path(regulation_dir, paste0("K", k, "_program", program_index, "_perturb_effects.txt"))
    tmp <- data.table::fread(program_path, data.table = FALSE)
    if (ncol(tmp) < 3) {
      stop(sprintf("Regulation file %s must have at least 3 columns", program_path))
    }
    colnames(tmp)[1:3] <- c("gene", "beta", "p")
    tmp$Program <- normalize_program_ids(program_index)
    rows[[length(rows) + 1]] <- tmp[, c("Program", "gene", "beta", "p")]
  }

  if (length(rows) == 0) {
    return(data.frame(
      Program = character(),
      gene = character(),
      beta = numeric(),
      p = numeric(),
      fdr = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  regulator_df <- do.call(rbind, rows)
  regulator_df$fdr <- stats::ave(
    regulator_df$p,
    regulator_df$Program,
    FUN = function(values) stats::p.adjust(values, method = "BH")
  )
  regulator_df
}

truncate_gene_label <- function(label, max_chars = 16) {
  label <- as.character(label)
  too_long <- nchar(label) > max_chars
  label[too_long] <- paste0(substr(label[too_long], 1, max_chars - 1), "\u2026")
  label
}
