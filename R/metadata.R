.read_metadata <- function(path) {
  # sample_id is always read as text so that IDs such as "001" keep their
  # leading zeros and still match RSEM file names.
  md <- .read_table_auto(path, character_cols = "sample_id")
  if (!"sample_id" %in% names(md)) {
    .stopf("Metadata must contain a column named 'sample_id'.")
  }
  md$sample_id <- as.character(md$sample_id)
  if (anyNA(md$sample_id) || any(!nzchar(md$sample_id))) {
    .stopf("Metadata contains missing or empty sample_id values.")
  }
  if (anyDuplicated(md$sample_id)) {
    dup <- unique(md$sample_id[duplicated(md$sample_id)])
    .stopf("Metadata contains duplicated sample_id values: %s", paste(dup, collapse = ", "))
  }
  md
}

.coerce_metadata <- function(md, categorical = character(), numeric = character()) {
  if (!is.data.frame(md)) .stopf("Metadata must be a data.frame.")
  missing_cat <- setdiff(categorical, names(md))
  missing_num <- setdiff(numeric, names(md))
  if (length(missing_cat)) .stopf("Unknown --p-categorical columns: %s", paste(missing_cat, collapse = ", "))
  if (length(missing_num)) .stopf("Unknown --p-numeric columns: %s", paste(missing_num, collapse = ", "))
  overlap <- intersect(categorical, numeric)
  if (length(overlap)) .stopf("Columns cannot be both categorical and numeric: %s", paste(overlap, collapse = ", "))

  out <- md
  for (nm in setdiff(names(out), "sample_id")) {
    x <- out[[nm]]
    if (nm %in% categorical) {
      out[[nm]] <- factor(as.character(x))
      next
    }
    if (nm %in% numeric) {
      y <- suppressWarnings(as.numeric(as.character(x)))
      if (any(is.na(y) & !is.na(x))) {
        .stopf("Column '%s' cannot be converted fully to numeric.", nm)
      }
      out[[nm]] <- y
      next
    }

    if (is.numeric(x) || is.integer(x)) {
      out[[nm]] <- as.numeric(x)
    } else if (is.logical(x)) {
      out[[nm]] <- factor(x)
    } else {
      chr <- as.character(x)
      y <- suppressWarnings(as.numeric(chr))
      non_missing <- !is.na(chr) & nzchar(chr)
      numeric_ok <- any(non_missing) && all(!is.na(y[non_missing]))
      if (numeric_ok) {
        out[[nm]] <- y
      } else {
        out[[nm]] <- factor(chr)
      }
    }
  }
  rownames(out) <- out$sample_id
  out
}

.apply_reference_levels <- function(md, reference = character()) {
  if (!length(reference)) return(md)
  out <- md
  for (spec in reference) {
    p <- .parse_double_colon(spec, 2L, "--p-reference")
    variable <- p[[1]]
    level <- p[[2]]
    if (!variable %in% names(out)) {
      .stopf("Reference variable '%s' is not present in metadata.", variable)
    }
    if (!is.factor(out[[variable]])) {
      out[[variable]] <- factor(out[[variable]])
    }
    if (!level %in% levels(out[[variable]])) {
      .stopf(
        "Reference level '%s' is not present in metadata column '%s'. Levels: %s",
        level, variable, paste(levels(out[[variable]]), collapse = ", ")
      )
    }
    out[[variable]] <- stats::relevel(out[[variable]], ref = level)
  }
  out
}

.formula_variables <- function(formula) {
  all.vars(formula)
}

.validate_formula_metadata <- function(md, formula) {
  vars <- .formula_variables(formula)
  missing <- setdiff(vars, names(md))
  if (length(missing)) {
    .stopf("Formula refers to metadata columns that do not exist: %s", paste(missing, collapse = ", "))
  }
  if (anyNA(md[, vars, drop = FALSE])) {
    bad <- vars[vapply(md[, vars, drop = FALSE], anyNA, logical(1))]
    .stopf(
      "Missing metadata values are not allowed in model variables: %s",
      paste(bad, collapse = ", ")
    )
  }
  invisible(TRUE)
}

.validate_design_matrix <- function(md, formula) {
  mm <- stats::model.matrix(formula, data = md)
  rank <- qr(mm)$rank
  if (rank < ncol(mm)) {
    .stopf(
      paste0(
        "The design matrix is not full rank (%d/%d independent columns). ",
        "The requested effects cannot be uniquely estimated. ",
        "Check for perfect confounding, redundant metadata columns, or empty factor combinations."
      ),
      rank, ncol(mm)
    )
  }
  mm
}
