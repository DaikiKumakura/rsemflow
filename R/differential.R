.prepare_deseq_metadata <- function(study, metadata = NULL, categorical = character(), numeric = character(), reference = character()) {
  md <- if (is.null(metadata)) {
    study$metadata
  } else if (is.character(metadata) && length(metadata) == 1L) {
    x <- .read_metadata(metadata)
    x <- x[match(study$sample_ids, x$sample_id), , drop = FALSE]
    if (anyNA(x$sample_id)) {
      .stopf("Metadata override does not contain all study samples.")
    }
    x
  } else if (is.data.frame(metadata)) {
    metadata
  } else {
    .stopf("metadata must be NULL, a metadata file path, or a data.frame.")
  }

  md <- .coerce_metadata(md, categorical = categorical, numeric = numeric)
  if (!all(study$sample_ids %in% md$sample_id)) {
    .stopf("Metadata is missing one or more samples from the study.")
  }
  md <- md[match(study$sample_ids, md$sample_id), , drop = FALSE]
  rownames(md) <- md$sample_id
  .apply_reference_levels(md, reference)
}

.find_numeric_coefficient <- function(dds, term) {
  rn <- DESeq2::resultsNames(dds)
  if (term %in% rn) return(term)
  exact <- rn[grepl(paste0("(^|_)", gsub("([.])", "\\\\\\1", term), "($|_)"), rn)]
  exact <- exact[!grepl("Intercept", exact, ignore.case = TRUE)]
  if (length(exact) == 1L) return(exact)
  .stopf(
    "Could not identify a unique DESeq2 coefficient for numeric term '%s'. Available coefficients: %s",
    term, paste(rn, collapse = ", ")
  )
}

.interaction_coefficients <- function(dds, term) {
  vars <- strsplit(term, ":", fixed = TRUE)[[1]]
  rn <- DESeq2::resultsNames(dds)
  keep <- vapply(rn, function(nm) all(vapply(vars, function(v) grepl(v, nm, fixed = TRUE), logical(1))), logical(1))
  hits <- rn[keep & !grepl("Intercept", rn, ignore.case = TRUE)]
  if (!length(hits)) {
    .stopf(
      "No DESeq2 interaction coefficients matched term '%s'. Available coefficients: %s",
      term, paste(rn, collapse = ", ")
    )
  }
  hits
}

.result_table <- function(res, gene_ids, term, contrast, coefficient = NA_character_, test = "Wald") {
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df$ensembl_gene_id <- .strip_ensembl_version(df$gene_id)
  df$term <- term
  df$contrast <- contrast
  df$coefficient <- coefficient
  df$test <- test
  front <- c("gene_id", "ensembl_gene_id", "term", "contrast", "coefficient", "test")
  df <- df[, c(front, setdiff(names(df), front)), drop = FALSE]
  rownames(df) <- NULL
  df
}

.parse_contrast <- function(spec) {
  p <- .parse_double_colon(spec, 3L, "--p-contrast")
  list(variable = p[[1]], numerator = p[[2]], denominator = p[[3]])
}

.simple_effect_at <- function(dds, md, contrast_spec, at_spec) {
  con <- .parse_contrast(contrast_spec)
  at <- .parse_double_colon(at_spec, 2L, "--p-at")
  at_var <- at[[1]]
  at_level <- at[[2]]
  var <- con$variable

  if (!is.factor(md[[var]]) || !is.factor(md[[at_var]])) {
    .stopf("--p-at currently requires categorical contrast and conditioning variables.")
  }
  ref_var <- levels(md[[var]])[[1]]
  ref_at <- levels(md[[at_var]])[[1]]
  if (con$denominator != ref_var) {
    .stopf(
      "--p-at currently requires the contrast denominator to be the reference level '%s' for '%s'.",
      ref_var, var
    )
  }
  if (!con$numerator %in% levels(md[[var]])) {
    .stopf("Unknown level '%s' for '%s'.", con$numerator, var)
  }
  if (!at_level %in% levels(md[[at_var]])) {
    .stopf("Unknown level '%s' for '%s'.", at_level, at_var)
  }

  main_res <- DESeq2::results(dds, contrast = c(var, con$numerator, con$denominator))
  if (at_level == ref_at) return(main_res)

  rn <- DESeq2::resultsNames(dds)
  candidates <- rn[
    grepl(var, rn, fixed = TRUE) &
      grepl(con$numerator, rn, fixed = TRUE) &
      grepl(at_var, rn, fixed = TRUE) &
      grepl(at_level, rn, fixed = TRUE)
  ]
  candidates <- candidates[!grepl("_vs_", candidates, fixed = TRUE)]
  if (length(candidates) != 1L) {
    .stopf(
      paste0(
        "Could not identify a unique interaction coefficient for '%s at %s::%s'. ",
        "Available coefficients: %s"
      ),
      contrast_spec, at_var, at_level, paste(rn, collapse = ", ")
    )
  }

  main_candidates <- rn[
    grepl(var, rn, fixed = TRUE) &
      grepl(con$numerator, rn, fixed = TRUE) &
      grepl("_vs_", rn, fixed = TRUE)
  ]
  if (length(main_candidates) != 1L) {
    .stopf(
      "Could not identify the main-effect coefficient for '%s'. Available coefficients: %s",
      contrast_spec, paste(rn, collapse = ", ")
    )
  }

  DESeq2::results(dds, contrast = list(c(main_candidates[[1]], candidates[[1]])))
}

#' Run metadata-driven DESeq2 differential expression
#'
#' @param study A `rsemflow_study` or study path.
#' @param formula R formula or formula string.
#' @param term Metadata variable or interaction term to report.
#' @param metadata Optional replacement metadata path/data.frame.
#' @param reference Character vector like `condition::Control`.
#' @param contrast Character vector like `condition::Drug::Control`.
#' @param at Optional conditioning specification like `genotype::KO`; currently
#'   supports one simple-effect condition with an explicit contrast.
#' @param test `wald` or `lrt`.
#' @param reduced Reduced formula required for LRT.
#' @param categorical Metadata columns to force categorical.
#' @param numeric Metadata columns to force numeric.
#' @param alpha Significance threshold passed to DESeq2 results.
#' @return A list with DESeq2 object, design matrix, contrast manifest, and
#'   named result tables.
#' @export
differential_deseq2 <- function(
    study,
    formula,
    term,
    metadata = NULL,
    reference = character(),
    contrast = character(),
    at = character(),
    test = c("wald", "lrt"),
    reduced = NULL,
    categorical = character(),
    numeric = character(),
    alpha = 0.05) {

  if (is.character(study) && length(study) == 1L) study <- read_study(study)
  test <- match.arg(tolower(test), c("wald", "lrt"))
  if (missing(formula) || is.null(formula)) .stopf("formula is required.")
  if (missing(term) || !nzchar(term)) .stopf("term is required.")

  f <- if (inherits(formula, "formula")) formula else stats::as.formula(as.character(formula))
  md <- .prepare_deseq_metadata(
    study, metadata = metadata, categorical = categorical,
    numeric = numeric, reference = reference
  )
  .validate_formula_metadata(md, f)
  design_matrix <- .validate_design_matrix(md, f)

  term_vars <- strsplit(term, ":", fixed = TRUE)[[1]]
  missing_term <- setdiff(term_vars, names(md))
  if (length(missing_term)) {
    .stopf("Term refers to unknown metadata columns: %s", paste(missing_term, collapse = ", "))
  }

  dds <- DESeq2::DESeqDataSetFromTximport(study$txi, colData = md, design = f)

  tables <- list()
  manifest <- data.frame(
    name = character(),
    term = character(),
    contrast = character(),
    coefficient = character(),
    stringsAsFactors = FALSE
  )

  if (test == "lrt") {
    if (is.null(reduced) || !length(reduced)) {
      .stopf("--p-reduced is required when --p-test lrt is used.")
    }
    rf <- if (inherits(reduced, "formula")) reduced else stats::as.formula(as.character(reduced))
    .validate_formula_metadata(md, rf)
    .validate_design_matrix(md, rf)

    dds <- DESeq2::DESeq(dds, test = "LRT", reduced = rf, quiet = TRUE)
    res <- DESeq2::results(dds, alpha = alpha)
    nm <- paste0(.safe_name(term), "_omnibus")
    lrt_tab <- .result_table(res, rownames(dds), term, "omnibus", NA_character_, test = "LRT")
    if ("log2FoldChange" %in% names(lrt_tab)) lrt_tab$log2FoldChange <- NA_real_
    if ("lfcSE" %in% names(lrt_tab)) lrt_tab$lfcSE <- NA_real_
    tables[[nm]] <- lrt_tab
    manifest <- rbind(
      manifest,
      data.frame(name = nm, term = term, contrast = "omnibus", coefficient = NA_character_, stringsAsFactors = FALSE)
    )
    return(list(
      dds = dds,
      metadata = md,
      formula = f,
      reduced = rf,
      test = test,
      design_matrix = design_matrix,
      manifest = manifest,
      results = tables
    ))
  }

  dds <- DESeq2::DESeq(dds, test = "Wald", quiet = TRUE)

  # Optional conditional simple effect. Kept intentionally narrow and explicit.
  if (length(at)) {
    if (length(at) != 1L || length(contrast) != 1L) {
      .stopf("--p-at currently requires exactly one --p-at and one --p-contrast.")
    }
    con <- .parse_contrast(contrast[[1]])
    atp <- .parse_double_colon(at[[1]], 2L, "--p-at")
    res <- .simple_effect_at(dds, md, contrast[[1]], at[[1]])
    label <- sprintf(
      "%s_%s_vs_%s_at_%s_%s",
      con$variable, con$numerator, con$denominator, atp[[1]], atp[[2]]
    )
    nm <- .safe_name(label)
    tables[[nm]] <- .result_table(
      res, rownames(dds), term,
      sprintf("%s vs %s at %s=%s", con$numerator, con$denominator, atp[[1]], atp[[2]]),
      NA_character_
    )
    manifest <- rbind(
      manifest,
      data.frame(
        name = nm, term = term,
        contrast = sprintf("%s::%s::%s @ %s::%s", con$variable, con$numerator, con$denominator, atp[[1]], atp[[2]]),
        coefficient = NA_character_,
        stringsAsFactors = FALSE
      )
    )
    return(list(
      dds = dds,
      metadata = md,
      formula = f,
      reduced = NULL,
      test = test,
      design_matrix = design_matrix,
      manifest = manifest,
      results = tables
    ))
  }

  if (grepl(":", term, fixed = TRUE)) {
    coeffs <- .interaction_coefficients(dds, term)
    for (coef in coeffs) {
      res <- DESeq2::results(dds, name = coef, alpha = alpha)
      nm <- .safe_name(coef)
      tables[[nm]] <- .result_table(res, rownames(dds), term, coef, coef)
      manifest <- rbind(
        manifest,
        data.frame(name = nm, term = term, contrast = coef, coefficient = coef, stringsAsFactors = FALSE)
      )
    }
  } else if (is.factor(md[[term]])) {
    lv <- levels(md[[term]])
    ref <- lv[[1]]
    for (level in setdiff(lv, ref)) {
      res <- DESeq2::results(dds, contrast = c(term, level, ref), alpha = alpha)
      label <- sprintf("%s_%s_vs_%s", term, level, ref)
      nm <- .safe_name(label)
      tables[[nm]] <- .result_table(
        res, rownames(dds), term, sprintf("%s vs %s", level, ref), NA_character_
      )
      manifest <- rbind(
        manifest,
        data.frame(name = nm, term = term, contrast = sprintf("%s::%s::%s", term, level, ref), coefficient = NA_character_, stringsAsFactors = FALSE)
      )
    }
  } else {
    coef <- .find_numeric_coefficient(dds, term)
    res <- DESeq2::results(dds, name = coef, alpha = alpha)
    nm <- .safe_name(term)
    tables[[nm]] <- .result_table(res, rownames(dds), term, "per_unit", coef)
    manifest <- rbind(
      manifest,
      data.frame(name = nm, term = term, contrast = "per_unit", coefficient = coef, stringsAsFactors = FALSE)
    )
  }

  if (length(contrast)) {
    for (spec in contrast) {
      con <- .parse_contrast(spec)
      if (!con$variable %in% names(md)) {
        .stopf("Contrast variable '%s' is not in metadata.", con$variable)
      }
      if (!is.factor(md[[con$variable]])) {
        .stopf("Explicit contrasts require a categorical metadata variable: %s", con$variable)
      }
      res <- DESeq2::results(
        dds,
        contrast = c(con$variable, con$numerator, con$denominator),
        alpha = alpha
      )
      label <- sprintf("%s_%s_vs_%s", con$variable, con$numerator, con$denominator)
      nm <- .safe_name(label)
      tables[[nm]] <- .result_table(
        res, rownames(dds), con$variable,
        sprintf("%s vs %s", con$numerator, con$denominator),
        NA_character_
      )
      manifest <- rbind(
        manifest,
        data.frame(name = nm, term = con$variable, contrast = spec, coefficient = NA_character_, stringsAsFactors = FALSE)
      )
    }
  }

  list(
    dds = dds,
    metadata = md,
    formula = f,
    reduced = NULL,
    test = test,
    design_matrix = design_matrix,
    manifest = unique(manifest),
    results = tables
  )
}
