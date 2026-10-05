.cli_actions <- function() {
  list(
    "data import" = list(
      summary = "Import RSEM *.genes.results files and sample metadata.",
      options = list(
        "i-rsem-dir" = c("required", "Directory containing *.genes.results files."),
        "m-metadata" = c("required", "TSV/CSV metadata containing sample_id."),
        "p-categorical" = c("multiple", "Metadata column to force categorical; repeatable."),
        "p-numeric" = c("multiple", "Metadata column to force numeric; repeatable."),
        "o-study" = c("required", "Output study directory.")
      )
    ),
    "data summarize" = list(
      summary = "Compute sample- and gene-level expression summaries.",
      options = list(
        "i-study" = c("required", "Input study directory or study.rds."),
        "o-summary" = c("required", "Output directory.")
      )
    ),
    "annotation ensembl" = list(
      summary = "Map versioned Ensembl gene IDs to stable IDs and gene symbols using a local OrgDb.",
      options = list(
        "i-study" = c("required", "Input study directory or study.rds."),
        "p-species" = c("default:auto", "auto, human, mouse, or rat."),
        "o-annotation" = c("required", "Output annotation TSV.")
      )
    ),
    "annotation gtf" = list(
      summary = "Annotate genes from a GTF/GTF.GZ file.",
      options = list(
        "i-study" = c("required", "Input study directory or study.rds."),
        "i-gtf" = c("required", "Input GTF or GTF.GZ."),
        "o-annotation" = c("required", "Output annotation TSV.")
      )
    ),
    "expression normalize" = list(
      summary = "Write a normalized/transformed expression matrix.",
      options = list(
        "i-study" = c("required", "Input study directory or study.rds."),
        "p-method" = c("default:vst", "vst, normalized-counts, or log2-tpm."),
        "o-expression" = c("required", "Output TSV.")
      )
    ),
    "expression pca" = list(
      summary = "Compute PCA scores, loadings, and explained variance.",
      options = list(
        "i-study" = c("required", "Input study directory or study.rds."),
        "p-transform" = c("default:vst", "vst or log2-tpm."),
        "p-ntop" = c("default:500", "Number of most variable genes."),
        "p-components" = c("default:10", "Maximum PCs to write."),
        "o-pca" = c("required", "Output directory.")
      )
    ),
    "expression correlation" = list(
      summary = "Compute a sample-by-sample expression correlation matrix.",
      options = list(
        "i-study" = c("required", "Input study directory or study.rds."),
        "p-transform" = c("default:vst", "vst or log2-tpm."),
        "p-method" = c("default:pearson", "pearson or spearman."),
        "o-correlation" = c("required", "Output TSV.")
      )
    ),
    "differential deseq2" = list(
      summary = "Run metadata-driven DESeq2 differential expression.",
      options = list(
        "i-study" = c("required", "Input study directory or study.rds."),
        "m-metadata" = c("optional", "Optional replacement TSV/CSV metadata."),
        "p-formula" = c("required", "Design formula, e.g. '~ batch + condition'."),
        "p-term" = c("required", "Metadata term to report/test."),
        "p-reference" = c("multiple", "Reference as COLUMN::LEVEL; repeatable."),
        "p-contrast" = c("multiple", "Explicit comparison COLUMN::NUMERATOR::DENOMINATOR; repeatable."),
        "p-at" = c("multiple", "Condition a simple effect as COLUMN::LEVEL; currently one value."),
        "p-categorical" = c("multiple", "Metadata column to force categorical; repeatable."),
        "p-numeric" = c("multiple", "Metadata column to force numeric; repeatable."),
        "p-test" = c("default:wald", "wald or lrt."),
        "p-reduced" = c("optional", "Reduced formula required for lrt."),
        "p-alpha" = c("default:0.05", "DESeq2 alpha threshold."),
        "o-differential" = c("required", "Output directory.")
      )
    ),
    "enrichment gsea" = list(
      summary = "Run preranked GSEA using a differential-expression statistic.",
      options = list(
        "i-differential" = c("required", "Input differential-expression TSV."),
        "i-genesets" = c("optional", "Optional custom GMT file."),
        "p-species" = c("default:auto", "auto, human, mouse, or rat."),
        "p-collection" = c("default:hallmark", "hallmark, reactome, go-bp, go-mf, go-cc, or MSigDB code."),
        "p-subcollection" = c("optional", "Raw MSigDB subcollection code."),
        "p-rank" = c("default:stat", "Differential table column used for ranking."),
        "p-min-size" = c("default:15", "Minimum mapped gene-set size."),
        "p-max-size" = c("default:500", "Maximum mapped gene-set size."),
        "p-eps" = c("default:0", "fgseaMultilevel eps."),
        "p-seed" = c("default:42", "Random seed for fgseaMultilevel."),
        "o-enrichment" = c("required", "Output TSV.")
      )
    ),
    "enrichment gsva" = list(
      summary = "Compute sample-level GSVA or ssGSEA pathway scores from log2(TPM+1).",
      options = list(
        "i-study" = c("required", "Input study directory or study.rds."),
        "i-genesets" = c("optional", "Optional custom GMT file."),
        "p-species" = c("default:auto", "auto, human, mouse, or rat."),
        "p-collection" = c("default:hallmark", "hallmark, reactome, go-bp, go-mf, go-cc, or MSigDB code."),
        "p-subcollection" = c("optional", "Raw MSigDB subcollection code."),
        "p-method" = c("default:gsva", "gsva or ssgsea."),
        "p-min-size" = c("default:10", "Minimum mapped gene-set size."),
        "p-max-size" = c("default:500", "Maximum mapped gene-set size."),
        "o-scores" = c("required", "Output TSV.")
      )
    )
  )
}

.parse_cli_options <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (token %in% c("--help", "-h")) {
      out$help <- TRUE
      i <- i + 1L
      next
    }
    if (!startsWith(token, "--")) {
      .stopf("Unexpected positional argument: %s", token)
    }
    key <- sub("^--", "", token)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      .stopf("Option --%s requires a value.", key)
    }
    value <- args[[i + 1L]]
    i <- i + 2L
    if (is.null(out[[key]])) out[[key]] <- value else out[[key]] <- c(out[[key]], value)
  }
  out
}

.opt <- function(opts, key, default = NULL) {
  if (is.null(opts[[key]])) default else opts[[key]]
}

.validate_cli_options <- function(action, opts) {
  spec <- .cli_actions()[[action]]
  if (is.null(spec)) .stopf("Unknown action: %s", action)
  allowed <- names(spec$options)
  unknown <- setdiff(names(opts), c(allowed, "help"))
  if (length(unknown)) {
    .stopf("Unknown option(s) for '%s': --%s", action, paste(unknown, collapse = ", --"))
  }
  for (nm in allowed) {
    type <- spec$options[[nm]][[1]]
    if (identical(type, "required") && is.null(opts[[nm]])) {
      .stopf("Missing required option --%s for '%s'.", nm, action)
    }
    if (!identical(type, "multiple") && length(opts[[nm]]) > 1L) {
      .stopf("Option --%s may be given only once.", nm)
    }
  }
  invisible(TRUE)
}

.opt_number <- function(opts, key, integer = FALSE) {
  raw <- .opt(opts, key)
  x <- suppressWarnings(as.numeric(raw))
  if (length(x) != 1L || is.na(x) || (integer && x != round(x))) {
    .stopf("Option --%s expects %s, got '%s'.", key, if (integer) "an integer" else "a number", raw)
  }
  if (integer) as.integer(x) else x
}

.option_default <- function(spec_entry) {
  tag <- spec_entry[[1]]
  if (!startsWith(tag, "default:")) return(NULL)
  sub("^default:", "", tag)
}

.print_main_help <- function() {
  cat(
"RSEMflow - metadata-driven downstream analysis for RSEM gene-level results

Usage:
  rsemflow <module> <action> [options]

Modules and actions:
  data          import, summarize
  annotation    ensembl, gtf
  expression    normalize, pca, correlation
  differential  deseq2
  enrichment    gsea, gsva

Option grammar:
  --i-*   input data
  --m-*   metadata
  --p-*   analysis parameters
  --o-*   outputs

Examples:
  rsemflow data import --help
  rsemflow differential deseq2 --help
  rsemflow enrichment gsea --help

Version:
  ", .rsemflow_version(), "\n", sep = "")
}

.print_module_help <- function(module) {
  actions <- .cli_actions()
  keys <- names(actions)[startsWith(names(actions), paste0(module, " "))]
  if (!length(keys)) .stopf("Unknown module: %s", module)
  cat(sprintf("RSEMflow module: %s\n\n", module))
  cat("Actions:\n")
  for (key in keys) {
    action <- sub(paste0("^", module, " "), "", key)
    cat(sprintf("  %-16s %s\n", action, actions[[key]]$summary))
  }
}

.print_action_help <- function(action) {
  spec <- .cli_actions()[[action]]
  if (is.null(spec)) .stopf("Unknown action: %s", action)
  parts <- strsplit(action, " ", fixed = TRUE)[[1]]
  cat(sprintf("%s\n\nUsage:\n  rsemflow %s %s [options]\n\n",
              spec$summary, parts[[1]], parts[[2]]))

  groups <- list(
    Inputs = grep("^i-", names(spec$options), value = TRUE),
    Metadata = grep("^m-", names(spec$options), value = TRUE),
    Parameters = grep("^p-", names(spec$options), value = TRUE),
    Outputs = grep("^o-", names(spec$options), value = TRUE)
  )
  for (grp in names(groups)) {
    opts <- groups[[grp]]
    if (!length(opts)) next
    cat(grp, ":\n", sep = "")
    for (nm in opts) {
      entry <- spec$options[[nm]]
      tag <- entry[[1]]
      desc <- entry[[2]]
      default <- .option_default(entry)
      suffix <- if (!is.null(default)) {
        sprintf(" [default: %s]", default)
      } else if (tag == "required") {
        " [required]"
      } else if (tag == "multiple") {
        " [repeatable]"
      } else ""
      cat(sprintf("  --%-22s %s%s\n", nm, desc, suffix))
    }
    cat("\n")
  }
}

.with_defaults <- function(action, opts) {
  spec <- .cli_actions()[[action]]
  for (nm in names(spec$options)) {
    if (is.null(opts[[nm]])) {
      default <- .option_default(spec$options[[nm]])
      if (!is.null(default)) opts[[nm]] <- default
    }
  }
  opts
}

.cli_data_import <- function(opts) {
  study <- read_rsem_study(
    rsem_dir = .opt(opts, "i-rsem-dir"),
    metadata = .opt(opts, "m-metadata"),
    categorical = .opt(opts, "p-categorical", character()),
    numeric = .opt(opts, "p-numeric", character())
  )
  write_study(study, .opt(opts, "o-study"))
  cat(sprintf(
    "Imported %d samples and %d genes into %s\n",
    length(study$sample_ids), nrow(study$txi$counts), .opt(opts, "o-study")
  ))
}

.cli_data_summarize <- function(opts) {
  out <- .opt(opts, "o-summary")
  s <- summarize_study(.opt(opts, "i-study"))
  .write_tsv(s$samples, file.path(out, "samples.tsv"))
  .write_tsv(s$genes, file.path(out, "genes.tsv"))
  cat(sprintf("Wrote expression summaries to %s\n", out))
}

.cli_annotation_ensembl <- function(opts) {
  ann <- annotate_ensembl(
    .opt(opts, "i-study"),
    species = .opt(opts, "p-species", "auto")
  )
  .write_tsv(ann, .opt(opts, "o-annotation"))
  cat(sprintf("Wrote %d gene annotations to %s\n", nrow(ann), .opt(opts, "o-annotation")))
}

.cli_annotation_gtf <- function(opts) {
  ann <- annotate_gtf(.opt(opts, "i-study"), .opt(opts, "i-gtf"))
  .write_tsv(ann, .opt(opts, "o-annotation"))
  cat(sprintf("Wrote %d gene annotations to %s\n", nrow(ann), .opt(opts, "o-annotation")))
}

.cli_expression_normalize <- function(opts) {
  mat <- normalize_expression(
    .opt(opts, "i-study"),
    method = .opt(opts, "p-method", "vst")
  )
  .write_tsv(.matrix_to_table(mat), .opt(opts, "o-expression"))
  cat(sprintf("Wrote expression matrix to %s\n", .opt(opts, "o-expression")))
}

.cli_expression_pca <- function(opts) {
  out <- .opt(opts, "o-pca")
  ntop <- .opt_number(opts, "p-ntop", integer = TRUE)
  components <- .opt_number(opts, "p-components", integer = TRUE)
  p <- compute_pca(
    .opt(opts, "i-study"),
    transform = .opt(opts, "p-transform", "vst"),
    ntop = ntop,
    components = components
  )
  .write_tsv(p$scores, file.path(out, "scores.tsv"))
  .write_tsv(p$loadings, file.path(out, "loadings.tsv"))
  .write_tsv(p$variance, file.path(out, "variance.tsv"))
  .write_key_value(
    list(transform = p$transform, ntop = p$ntop),
    file.path(out, "analysis-info.tsv")
  )
  cat(sprintf("Wrote PCA tables to %s\n", out))
}

.cli_expression_correlation <- function(opts) {
  mat <- compute_correlation(
    .opt(opts, "i-study"),
    transform = .opt(opts, "p-transform", "vst"),
    method = .opt(opts, "p-method", "pearson")
  )
  .write_tsv(.matrix_to_table(mat, id_name = "sample_id"), .opt(opts, "o-correlation"))
  cat(sprintf("Wrote sample correlation matrix to %s\n", .opt(opts, "o-correlation")))
}

.cli_differential_deseq2 <- function(opts) {
  out <- .opt(opts, "o-differential")
  alpha <- .opt_number(opts, "p-alpha")
  x <- differential_deseq2(
    study = .opt(opts, "i-study"),
    metadata = .opt(opts, "m-metadata"),
    formula = .opt(opts, "p-formula"),
    term = .opt(opts, "p-term"),
    reference = .opt(opts, "p-reference", character()),
    contrast = .opt(opts, "p-contrast", character()),
    at = .opt(opts, "p-at", character()),
    categorical = .opt(opts, "p-categorical", character()),
    numeric = .opt(opts, "p-numeric", character()),
    test = .opt(opts, "p-test", "wald"),
    reduced = .opt(opts, "p-reduced"),
    alpha = alpha
  )

  design <- data.frame(
    sample_id = rownames(x$design_matrix),
    as.data.frame(x$design_matrix, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  .write_tsv(design, file.path(out, "design-matrix.tsv"))
  .write_tsv(x$manifest, file.path(out, "contrasts.tsv"))
  for (nm in names(x$results)) {
    .write_tsv(x$results[[nm]], file.path(out, paste0(nm, ".tsv")))
  }
  info <- list(
    formula = paste(deparse(x$formula), collapse = ""),
    term = .opt(opts, "p-term"),
    test = x$test,
    reduced = if (is.null(x$reduced)) "" else paste(deparse(x$reduced), collapse = ""),
    alpha = x$alpha
  )
  .write_key_value(info, file.path(out, "analysis-info.tsv"))
  cat(sprintf("Wrote %d differential result table(s) to %s\n", length(x$results), out))
}

.cli_enrichment_gsea <- function(opts) {
  sizes <- c(.opt_number(opts, "p-min-size", integer = TRUE), .opt_number(opts, "p-max-size", integer = TRUE))
  eps <- .opt_number(opts, "p-eps")
  seed <- .opt_number(opts, "p-seed", integer = TRUE)
  res <- run_gsea(
    differential = .opt(opts, "i-differential"),
    genesets = .opt(opts, "i-genesets"),
    species = .opt(opts, "p-species", "auto"),
    collection = .opt(opts, "p-collection", "hallmark"),
    subcollection = .opt(opts, "p-subcollection"),
    rank = .opt(opts, "p-rank", "stat"),
    min_size = sizes[[1]],
    max_size = sizes[[2]],
    eps = eps,
    seed = seed
  )
  .write_tsv(res, .opt(opts, "o-enrichment"))
  cat(sprintf("Wrote %d GSEA results to %s\n", nrow(res), .opt(opts, "o-enrichment")))
}

.cli_enrichment_gsva <- function(opts) {
  sizes <- c(.opt_number(opts, "p-min-size", integer = TRUE), .opt_number(opts, "p-max-size", integer = TRUE))
  res <- run_gsva(
    study = .opt(opts, "i-study"),
    genesets = .opt(opts, "i-genesets"),
    species = .opt(opts, "p-species", "auto"),
    collection = .opt(opts, "p-collection", "hallmark"),
    subcollection = .opt(opts, "p-subcollection"),
    method = .opt(opts, "p-method", "gsva"),
    min_size = sizes[[1]],
    max_size = sizes[[2]]
  )
  .write_tsv(res, .opt(opts, "o-scores"))
  cat(sprintf("Wrote %d pathway score rows to %s\n", nrow(res), .opt(opts, "o-scores")))
}

.dispatch_cli <- function(action, opts) {
  switch(
    action,
    "data import" = .cli_data_import(opts),
    "data summarize" = .cli_data_summarize(opts),
    "annotation ensembl" = .cli_annotation_ensembl(opts),
    "annotation gtf" = .cli_annotation_gtf(opts),
    "expression normalize" = .cli_expression_normalize(opts),
    "expression pca" = .cli_expression_pca(opts),
    "expression correlation" = .cli_expression_correlation(opts),
    "differential deseq2" = .cli_differential_deseq2(opts),
    "enrichment gsea" = .cli_enrichment_gsea(opts),
    "enrichment gsva" = .cli_enrichment_gsva(opts),
    .stopf("Unknown action: %s", action)
  )
}

#' RSEMflow command-line entry point
#'
#' Parses `rsemflow <module> <action> [options]` and calls the matching R
#' function. Options use the prefixes `--i-` (inputs), `--m-` (metadata),
#' `--p-` (parameters), and `--o-` (outputs). Add `--help` at any level for
#' usage.
#'
#' The installed launcher script (`system.file("exec", "rsemflow", package =
#' "rsemflow")`) calls this function. Without the launcher, use
#' `Rscript -e "rsemflow::cli_main()" <module> <action> [options]`.
#'
#' @param args Character vector of command-line arguments.
#' @return Invisibly returns zero on success. Errors are raised as R errors.
#' @examples
#' cli_main("--help")
#' cli_main(c("differential", "deseq2", "--help"))
#' @export
cli_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || args[[1]] %in% c("--help", "-h")) {
    .print_main_help()
    return(invisible(0L))
  }
  if (args[[1]] %in% c("--version", "-V")) {
    cat(.rsemflow_version(), "\n", sep = "")
    return(invisible(0L))
  }

  module <- args[[1]]
  if (length(args) == 1L || args[[2]] %in% c("--help", "-h")) {
    .print_module_help(module)
    return(invisible(0L))
  }

  action_name <- paste(module, args[[2]])
  if (is.null(.cli_actions()[[action_name]])) {
    .stopf("Unknown module/action: %s", action_name)
  }

  opts <- .parse_cli_options(args[-c(1, 2)])
  if (isTRUE(opts$help)) {
    .print_action_help(action_name)
    return(invisible(0L))
  }
  .validate_cli_options(action_name, opts)
  opts <- .with_defaults(action_name, opts)
  .dispatch_cli(action_name, opts)
  invisible(0L)
}
