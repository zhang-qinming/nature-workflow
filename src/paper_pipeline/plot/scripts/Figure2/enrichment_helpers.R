has_rows <- function(x) {
  !is.null(x) && nrow(as.data.frame(x)) > 0
}

write_enrichment_table <- function(summary_df, output_path) {
  if (is.null(summary_df) || nrow(summary_df) == 0) {
    summary_df <- data.frame(
      TOP_N = character(),
      Description = character(),
      pvalue = numeric(),
      stringsAsFactors = FALSE
    )
  }

  if (!"TOP_N" %in% colnames(summary_df)) {
    summary_df$TOP_N <- character(nrow(summary_df))
  }
  if (!"Description" %in% colnames(summary_df)) {
    summary_df$Description <- character(nrow(summary_df))
  }
  if (!"pvalue" %in% colnames(summary_df)) {
    summary_df$pvalue <- numeric(nrow(summary_df))
  }

  write.table(
    data.frame(ID = row.names(summary_df), summary_df, row.names = NULL, check.names = FALSE),
    output_path,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}

load_hallmark_term2gene <- function(geneset_dir = file.path("data", "geneset")) {
  hallmark_files <- sort(Sys.glob(file.path(geneset_dir, "HALLMARK_*.txt")))

  if (length(hallmark_files) >= 50) {
    message(
      "Using ",
      length(hallmark_files),
      " local Hallmark gene-set files from `",
      geneset_dir,
      "`."
    )

    local_sets <- lapply(hallmark_files, function(path) {
      genes <- read.table(
        path,
        header = FALSE,
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = ""
      )[, 1]
      genes <- unique(stats::na.omit(genes))
      if (length(genes) == 0) {
        return(NULL)
      }

      mapped <- clusterProfiler::bitr(
        genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db,
        drop = TRUE
      )
      if (!has_rows(mapped)) {
        return(NULL)
      }

      mapped <- mapped[!duplicated(mapped$ENTREZID), "ENTREZID", drop = FALSE]
      data.frame(
        gs_name = tools::file_path_sans_ext(basename(path)),
        entrez_gene = mapped$ENTREZID,
        stringsAsFactors = FALSE
      )
    })

    local_sets <- Filter(Negate(is.null), local_sets)
    if (length(local_sets) > 0) {
      return(as.data.frame(unique(data.table::rbindlist(local_sets, use.names = TRUE, fill = TRUE))))
    }
  }

  if (length(hallmark_files) > 0) {
    warning(
      "Found only ",
      length(hallmark_files),
      " local `HALLMARK_*.txt` files under `",
      geneset_dir,
      "`. A full offline Hallmark collection should include 50 files, so the script will fall back to `msigdbr`."
    )
  }

  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop(
      "Failed to load Hallmark gene sets: no usable local `HALLMARK_*.txt` files were found under `",
      geneset_dir,
      "`, and the `msigdbr` package is not available.",
      call. = FALSE
    )
  }

  msigdbr_fun <- get("msigdbr", asNamespace("msigdbr"))
  msigdbr_args <- names(formals(msigdbr_fun))
  query <- list(species = "Homo sapiens")
  if ("collection" %in% msigdbr_args) {
    query$collection <- "H"
  } else {
    query$category <- "H"
  }

  tryCatch(
    {
      msig <- do.call(msigdbr_fun, query)
      as.data.frame(unique(msig[, c("gs_name", "entrez_gene")]))
    },
    error = function(e) {
      stop(
        "Failed to load Hallmark gene sets. Provide local `HALLMARK_*.txt` files under `",
        geneset_dir,
        "` via `geneset_dir` (50 files for the full Hallmark collection), or run `msigdbr` with a populated cache on the cluster. Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}
