.prepare_deseq_metadata <- function(study, metadata = NULL, categorical = character(), numeric = character(), reference = character()) {
  md <- if (is.null(metadata)) {
    study$metadata
  } else if (is.character(metadata) && length(metadata) == 1L) {
    .read_metadata(metadata)
  } else if (is.data.frame(metadata)) {
    metadata
  } else {
    .stopf("metadata must be NULL, a metadata file path, or a data.frame.")
  }
  if (!"sample_id" %in% names(md)) {
    .stopf("Metadata must contain a column named 'sample_id'.")
  }
  md$sample_id <- as.character(md$sample_id)
  missing <- setdiff(study$sample_ids, md$sample_id)
  if (length(missing)) {
    .stopf("Metadata is missing study samples: %s", paste(missing, collapse = ", "))
  }

  md <- md[match(study$sample_ids, md$sample_id), , drop = FALSE]
  md <- .coerce_metadata(md, categorical = categorical, numeric = numeric)
  rownames(md) <- md$sample_id
  .apply_reference_levels(md, reference)
}

.sorted_vars <- function(label) {
  paste(sort(strsplit(label, ":", fixed = TRUE)[[1]]), collapse = ":")
}

# Model-matrix columns that belong to one formula term, named the way DESeq2
# names its coefficients (make.names() of the model-matrix column names).
.term_columns <- function(design_matrix, formula, term) {
  labels <- attr(stats::terms(formula), "term.labels")
  idx <- match(.sorted_vars(term), vapply(labels, .sorted_vars, character(1)))
  if (is.na(idx)) {
    .stopf(
      "Term '%s' is not part of the formula. Formula terms: %s",
      term, paste(labels, collapse = ", ")
    )
  }
  cols <- colnames(design_matrix)[attr(design_matrix, "assign") == idx]
  make.names(cols)
}

.require_coefficients <- function(dds, coefs, what) {
  rn <- DESeq2::resultsNames(dds)
  missing <- setdiff(coefs, rn)
  if (length(missing)) {
    .stopf(
      "Could not find DESeq2 coefficient(s) %s for %s. Available coefficients: %s",
      paste(missing, collapse = ", "), what, paste(rn, collapse = ", ")
    )
  }
  coefs
}

.result_table <- function(res, term, contrast, coefficient = NA_character_, test = "Wald") {
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

.check_contrast_levels <- function(md, con) {
  if (!con$variable %in% names(md)) {
    .stopf("Contrast variable '%s' is not in metadata.", con$variable)
  }
  if (!is.factor(md[[con$variable]])) {
    .stopf("Explicit contrasts require a categorical metadata variable: %s", con$variable)
  }
  lv <- levels(md[[con$variable]])
  bad <- setdiff(c(con$numerator, con$denominator), lv)
  if (length(bad)) {
    .stopf(
      "Unknown level(s) %s for '%s'. Levels: %s",
      paste(bad, collapse = ", "), con$variable, paste(lv, collapse = ", ")
    )
  }
  if (identical(con$numerator, con$denominator)) {
    .stopf("Contrast numerator and denominator must differ: %s", con$variable)
  }
  invisible(TRUE)
}

# Simple effect of `var` (numerator vs reference) at one level of `at_var`:
# main-effect coefficient + the matching interaction coefficient.
.simple_effect_at <- function(dds, md, design_matrix, formula, contrast_spec, at_spec, alpha) {
  con <- .parse_contrast(contrast_spec)
  at <- .parse_double_colon(at_spec, 2L, "--p-at")
  at_var <- at[[1]]
  at_level <- at[[2]]
  var <- con$variable

  .check_contrast_levels(md, con)
  if (!at_var %in% names(md) || !is.factor(md[[at_var]])) {
    .stopf("--p-at requires a categorical metadata variable: %s", at_var)
  }
  if (!at_level %in% levels(md[[at_var]])) {
    .stopf(
      "Unknown level '%s' for '%s'. Levels: %s",
      at_level, at_var, paste(levels(md[[at_var]]), collapse = ", ")
    )
  }
  ref_var <- levels(md[[var]])[[1]]
  if (con$denominator != ref_var) {
    .stopf(
      "--p-at currently requires the contrast denominator to be the reference level '%s' for '%s'.",
      ref_var, var
    )
  }

  main <- .require_coefficients(
    dds,
    make.names(sprintf("%s_%s_vs_%s", var, con$numerator, ref_var)),
    sprintf("the main effect of '%s'", var)
  )
  if (at_level == levels(md[[at_var]])[[1]]) {
    return(list(
      res = DESeq2::results(dds, name = main, alpha = alpha),
      coefficient = main
    ))
  }

  interaction <- paste(var, at_var, sep = ":")
  inter_cols <- .term_columns(design_matrix, formula, interaction)
  wanted <- make.names(c(
    paste0(var, con$numerator, ":", at_var, at_level),
    paste0(at_var, at_level, ":", var, con$numerator)
  ))
  inter <- intersect(wanted, inter_cols)
  if (length(inter) != 1L) {
    .stopf(
      "Could not identify the interaction coefficient for '%s' at %s::%s.",
      contrast_spec, at_var, at_level
    )
  }
  .require_coefficients(dds, inter, sprintf("the interaction '%s'", interaction))
  list(
    res = DESeq2::results(dds, contrast = list(c(main, inter)), alpha = alpha),
    coefficient = paste(main, inter, sep = "+")
  )
}

.interacting_terms <- function(formula, var) {
  labels <- attr(stats::terms(formula), "term.labels")
  labels[grepl(":", labels, fixed = TRUE) &
    vapply(strsplit(labels, ":", fixed = TRUE), function(v) var %in% v, logical(1))]
}

.manifest_row <- function(name, term, contrast, coefficient) {
  data.frame(
    name = name, term = term, contrast = contrast, coefficient = coefficient,
    stringsAsFactors = FALSE
  )
}

.differential_output <- function(dds, md, formula, reduced, test, alpha, design_matrix, manifest, tables) {
  list(
    dds = dds,
    metadata = md,
    formula = formula,
    reduced = reduced,
    test = test,
    alpha = alpha,
    design_matrix = design_matrix,
    manifest = manifest,
    results = tables
  )
}

#' Run metadata-driven DESeq2 differential expression
#'
#' Fits a DESeq2 model from the tximport-derived counts of a study and returns
#' one result table per comparison.
#'
#' Which tables are produced depends on `term`:
#'
#' * a categorical variable: every level versus the reference level;
#' * a numeric variable: the per-unit coefficient;
#' * an interaction such as `genotype:treatment`: every interaction coefficient.
#'
#' Extra comparisons can be requested with `contrast`. With `at`, a single
#' simple effect is returned instead (for example the treatment effect within
#' one genotype in a `~ genotype * treatment` model). With `test = "lrt"`, one
#' omnibus table is returned and `log2FoldChange`/`lfcSE` are set to `NA`
#' because the likelihood-ratio test has no single direction.
#'
#' Log fold changes are not shrunk.
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
#' @return A list with the fitted `DESeqDataSet` (`dds`), the metadata used,
#'   the formula(s), the design matrix, a contrast manifest, and the named
#'   result tables (`results`).
#' @examples
#' \donttest{
#' ex <- system.file("extdata", "example", package = "rsemflow")
#' study <- read_rsem_study(file.path(ex, "data"), file.path(ex, "metadata.tsv"))
#' de <- differential_deseq2(
#'   study,
#'   formula = ~ batch + condition,
#'   term = "condition",
#'   reference = "condition::Control"
#' )
#' head(de$results$condition_Drug_vs_Control)
#' }
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
  test <- match.arg(tolower(test[[1]]), c("wald", "lrt"))
  if (missing(formula) || is.null(formula)) .stopf("formula is required.")
  if (missing(term) || length(term) != 1L || !nzchar(term)) .stopf("term is required.")
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    .stopf("alpha must be a single number between 0 and 1.")
  }

  f <- if (inherits(formula, "formula")) formula else stats::as.formula(as.character(formula))
  md <- .prepare_deseq_metadata(
    study,
    metadata = metadata, categorical = categorical,
    numeric = numeric, reference = reference
  )
  .validate_formula_metadata(md, f)
  design_matrix <- .validate_design_matrix(md, f)

  term_vars <- strsplit(term, ":", fixed = TRUE)[[1]]
  missing_term <- setdiff(term_vars, names(md))
  if (length(missing_term)) {
    .stopf("Term refers to unknown metadata columns: %s", paste(missing_term, collapse = ", "))
  }
  missing_in_formula <- setdiff(term_vars, all.vars(f))
  if (length(missing_in_formula)) {
    .stopf("Term variable(s) not used in the formula: %s", paste(missing_in_formula, collapse = ", "))
  }

  dds <- .deseq_dataset(study, md, f)
  tables <- list()
  manifest <- .manifest_row(character(), character(), character(), character())

  if (test == "lrt") {
    if (is.null(reduced) || !length(reduced) || !nzchar(reduced[[1]])) {
      .stopf("--p-reduced is required when --p-test lrt is used.")
    }
    rf <- if (inherits(reduced, "formula")) reduced else stats::as.formula(as.character(reduced))
    .validate_formula_metadata(md, rf)
    .validate_design_matrix(md, rf)
    full_labels <- vapply(attr(stats::terms(f), "term.labels"), .sorted_vars, character(1))
    reduced_labels <- vapply(attr(stats::terms(rf), "term.labels"), .sorted_vars, character(1))
    if (length(setdiff(reduced_labels, full_labels))) {
      .stopf("The reduced formula must be nested in the full formula.")
    }
    dropped <- setdiff(full_labels, reduced_labels)
    if (!.sorted_vars(term) %in% dropped) {
      .stopf(
        "Term '%s' must be removed in the reduced formula. Terms removed: %s",
        term, if (length(dropped)) paste(dropped, collapse = ", ") else "(none)"
      )
    }

    dds <- .quietly(DESeq2::DESeq(dds, test = "LRT", reduced = rf, quiet = TRUE))
    res <- DESeq2::results(dds, alpha = alpha)
    nm <- paste0(.safe_name(term), "_omnibus")
    lrt_tab <- .result_table(res, term, "omnibus", NA_character_, test = "LRT")
    lrt_tab$log2FoldChange <- NA_real_
    lrt_tab$lfcSE <- NA_real_
    tables[[nm]] <- lrt_tab
    manifest <- .manifest_row(nm, term, "omnibus", NA_character_)
    return(.differential_output(dds, md, f, rf, test, alpha, design_matrix, manifest, tables))
  }

  dds <- .quietly(DESeq2::DESeq(dds, test = "Wald", quiet = TRUE))

  # Optional conditional simple effect. Kept intentionally narrow and explicit.
  if (length(at)) {
    if (length(at) != 1L || length(contrast) != 1L) {
      .stopf("--p-at currently requires exactly one --p-at and one --p-contrast.")
    }
    con <- .parse_contrast(contrast[[1]])
    atp <- .parse_double_colon(at[[1]], 2L, "--p-at")
    eff <- .simple_effect_at(dds, md, design_matrix, f, contrast[[1]], at[[1]], alpha)
    nm <- .safe_name(sprintf(
      "%s_%s_vs_%s_at_%s_%s",
      con$variable, con$numerator, con$denominator, atp[[1]], atp[[2]]
    ))
    tables[[nm]] <- .result_table(
      eff$res, term,
      sprintf("%s vs %s at %s=%s", con$numerator, con$denominator, atp[[1]], atp[[2]]),
      eff$coefficient
    )
    manifest <- .manifest_row(
      nm, term,
      sprintf("%s::%s::%s @ %s::%s", con$variable, con$numerator, con$denominator, atp[[1]], atp[[2]]),
      eff$coefficient
    )
    return(.differential_output(dds, md, f, NULL, test, alpha, design_matrix, manifest, tables))
  }

  if (length(term_vars) > 1L) {
    coeffs <- .require_coefficients(
      dds, .term_columns(design_matrix, f, term), sprintf("term '%s'", term)
    )
    for (coef in coeffs) {
      res <- DESeq2::results(dds, name = coef, alpha = alpha)
      nm <- .safe_name(coef)
      tables[[nm]] <- .result_table(res, term, coef, coef)
      manifest <- rbind(manifest, .manifest_row(nm, term, coef, coef))
    }
  } else if (is.factor(md[[term]])) {
    interacting <- .interacting_terms(f, term)
    if (length(interacting)) {
      others <- setdiff(unlist(strsplit(interacting, ":", fixed = TRUE)), term)
      .warnf(
        paste0(
          "'%s' interacts with %s in the formula, so the '%s' comparisons are ",
          "simple effects at the reference level(s) of %s. Use --p-at for other levels."
        ),
        term, paste(unique(others), collapse = ", "), term,
        paste(sprintf("%s (%s)", unique(others), vapply(unique(others), function(v) {
          if (is.factor(md[[v]])) levels(md[[v]])[[1]] else "0"
        }, character(1))), collapse = ", ")
      )
    }
    lv <- levels(md[[term]])
    ref <- lv[[1]]
    for (level in setdiff(lv, ref)) {
      res <- DESeq2::results(dds, contrast = c(term, level, ref), alpha = alpha)
      nm <- .safe_name(sprintf("%s_%s_vs_%s", term, level, ref))
      tables[[nm]] <- .result_table(res, term, sprintf("%s vs %s", level, ref), NA_character_)
      manifest <- rbind(
        manifest,
        .manifest_row(nm, term, sprintf("%s::%s::%s", term, level, ref), NA_character_)
      )
    }
  } else {
    coef <- .require_coefficients(
      dds, .term_columns(design_matrix, f, term), sprintf("numeric term '%s'", term)
    )
    res <- DESeq2::results(dds, name = coef, alpha = alpha)
    nm <- .safe_name(term)
    tables[[nm]] <- .result_table(res, term, "per_unit", coef)
    manifest <- rbind(manifest, .manifest_row(nm, term, "per_unit", coef))
  }

  for (spec in contrast) {
    con <- .parse_contrast(spec)
    .check_contrast_levels(md, con)
    nm <- .safe_name(sprintf("%s_%s_vs_%s", con$variable, con$numerator, con$denominator))
    if (nm %in% names(tables)) next
    res <- DESeq2::results(
      dds,
      contrast = c(con$variable, con$numerator, con$denominator),
      alpha = alpha
    )
    tables[[nm]] <- .result_table(
      res, con$variable,
      sprintf("%s vs %s", con$numerator, con$denominator),
      NA_character_
    )
    manifest <- rbind(manifest, .manifest_row(nm, con$variable, spec, NA_character_))
  }

  .differential_output(dds, md, f, NULL, test, alpha, design_matrix, manifest, tables)
}
