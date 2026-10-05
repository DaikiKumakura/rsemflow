.rsemflow_version <- function() {
  as.character(utils::packageVersion("rsemflow"))
}

.stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

.warnf <- function(fmt, ...) {
  warning(sprintf(fmt, ...), call. = FALSE)
}

# Bioconductor functions print progress messages ("reading in files",
# "using counts and average transcript lengths") that clutter CLI output.
# Warnings and errors are still propagated.
.quietly <- function(expr) {
  suppressMessages(expr)
}

.ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(path)) {
      .stopf("Could not create directory: %s", path)
    }
  }
  invisible(normalizePath(path, mustWork = FALSE))
}

.write_tsv <- function(x, path) {
  .ensure_dir(dirname(path))
  data.table::fwrite(
    as.data.frame(x, stringsAsFactors = FALSE),
    file = path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  invisible(path)
}

.read_table_auto <- function(path, character_cols = NULL) {
  if (!file.exists(path)) {
    .stopf("File does not exist: %s", path)
  }
  ext <- tolower(tools::file_ext(path))
  sep <- if (ext == "csv") "," else "\t"
  if (length(character_cols)) {
    header <- names(data.table::fread(path, sep = sep, nrows = 0L, check.names = FALSE))
    character_cols <- intersect(character_cols, header)
  }
  out <- data.table::fread(
    path,
    sep = sep,
    header = TRUE,
    colClasses = if (length(character_cols)) list(character = character_cols),
    encoding = "UTF-8",
    showProgress = FALSE,
    data.table = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "N/A", "NaN")
  )
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

.matrix_to_table <- function(mat, id_name = "gene_id") {
  if (is.null(rownames(mat))) {
    .stopf("Matrix has no row names and cannot be written as a feature table.")
  }
  out <- data.frame(
    .id = rownames(mat),
    as.data.frame(mat, check.names = FALSE, stringsAsFactors = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(out)[1] <- id_name
  out
}

# Only Ensembl gene IDs lose their version suffix. Other identifiers (GENCODE
# `_PAR_Y` IDs, or names such as `Y_RNA.1`) are kept unchanged so that
# unrelated rows are never merged.
.strip_ensembl_version <- function(x) {
  sub("^(ENS[A-Z]*G[0-9]+)\\.[0-9]+$", "\\1", as.character(x))
}

.looks_like_ensembl <- function(x) {
  x <- .strip_ensembl_version(x)
  grepl("^ENS[A-Z]*G[0-9]+$", x)
}

.infer_species_from_ids <- function(ids) {
  ids <- .strip_ensembl_version(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) return(NA_character_)
  props <- c(
    human = mean(grepl("^ENSG[0-9]+$", ids)),
    mouse = mean(grepl("^ENSMUSG[0-9]+$", ids)),
    rat = mean(grepl("^ENSRNOG[0-9]+$", ids))
  )
  best <- names(which.max(props))
  if (max(props) < 0.5) NA_character_ else best
}

.resolve_species <- function(species, ids = NULL) {
  species <- tolower(trimws(species %||% "auto"))
  if (species == "auto") {
    inferred <- .infer_species_from_ids(ids)
    if (is.na(inferred)) {
      .stopf("Species could not be inferred from gene IDs. Specify --p-species explicitly.")
    }
    species <- inferred
  }
  aliases <- list(
    human = c("human", "homo sapiens", "hs", "h.sapiens"),
    mouse = c("mouse", "mus musculus", "mm", "m.musculus"),
    rat = c("rat", "rattus norvegicus", "rn", "r.norvegicus")
  )
  for (nm in names(aliases)) {
    if (species %in% aliases[[nm]]) return(nm)
  }
  .stopf("Unsupported species '%s'. Supported values: auto, human, mouse, rat.", species)
}

.species_scientific <- function(species) {
  switch(
    .resolve_species(species),
    human = "Homo sapiens",
    mouse = "Mus musculus",
    rat = "Rattus norvegicus"
  )
}

.species_db <- function(species) {
  switch(
    .resolve_species(species),
    human = "HS",
    mouse = "MM",
    # MSigDB has no native rat database; msigdbr maps human sets to rat orthologs.
    rat = "HS"
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

.as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x)) return(default)
  if (is.logical(x)) return(x[[1]])
  y <- tolower(trimws(as.character(x[[1]])))
  if (y %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (y %in% c("false", "f", "0", "no", "n")) return(FALSE)
  .stopf("Expected a boolean value, got '%s'.", x[[1]])
}

.safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "result")
}

.parse_double_colon <- function(x, expected, option_name) {
  parts <- strsplit(as.character(x), "::", fixed = TRUE)[[1]]
  if (length(parts) != expected || any(!nzchar(parts))) {
    .stopf(
      "%s expects %d fields separated by '::'; got '%s'.",
      option_name, expected, x
    )
  }
  parts
}

.collapse_matrix_by_id <- function(mat, ids) {
  ids <- as.character(ids)
  keep <- !is.na(ids) & nzchar(ids)
  mat <- mat[keep, , drop = FALSE]
  ids <- ids[keep]
  if (!anyDuplicated(ids)) {
    rownames(mat) <- ids
    return(mat)
  }
  summed <- rowsum(mat, group = ids, reorder = FALSE, na.rm = TRUE)
  n <- as.numeric(table(factor(ids, levels = rownames(summed))))
  summed / n
}

.available_pc_count <- function(pca) {
  min(ncol(pca$x), ncol(pca$rotation))
}

.write_key_value <- function(x, path) {
  if (is.null(names(x))) .stopf("Key-value output requires a named vector/list.")
  vals <- vapply(x, function(v) paste(as.character(v), collapse = ","), character(1))
  .write_tsv(data.frame(key = names(vals), value = vals, stringsAsFactors = FALSE), path)
}
