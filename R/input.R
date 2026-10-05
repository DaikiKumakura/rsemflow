.find_rsem_files <- function(rsem_dir) {
  if (!dir.exists(rsem_dir)) {
    .stopf("RSEM input directory does not exist: %s", rsem_dir)
  }
  files <- list.files(
    rsem_dir,
    pattern = "\\.genes\\.results$",
    full.names = TRUE
  )
  if (!length(files)) {
    .stopf("No '*.genes.results' files were found in: %s", rsem_dir)
  }
  files <- sort(files)
  sample_ids <- sub("\\.genes\\.results$", "", basename(files))
  if (anyDuplicated(sample_ids)) {
    .stopf("Duplicate sample IDs were inferred from RSEM filenames.")
  }
  names(files) <- sample_ids
  files
}

.validate_rsem_header <- function(path) {
  x <- data.table::fread(
    path,
    nrows = 1L,
    data.table = FALSE,
    check.names = FALSE,
    showProgress = FALSE
  )
  required <- c(
    "gene_id", "transcript_id(s)", "length", "effective_length",
    "expected_count", "TPM", "FPKM"
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    .stopf(
      "RSEM file '%s' is missing required columns: %s",
      path, paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}

.read_gene_ids <- function(path) {
  x <- data.table::fread(
    path,
    select = "gene_id",
    data.table = FALSE,
    showProgress = FALSE
  )
  as.character(x$gene_id)
}

.validate_rsem_gene_ids <- function(files) {
  ref <- .read_gene_ids(files[[1]])
  if (anyDuplicated(ref)) {
    .stopf("RSEM file contains duplicated gene_id values: %s", files[[1]])
  }
  if (length(files) > 1L) {
    for (i in 2:length(files)) {
      ids <- .read_gene_ids(files[[i]])
      if (!identical(ids, ref)) {
        .stopf(
          paste0(
            "Gene IDs or their order differ between RSEM files '%s' and '%s'. ",
            "RSEMflow requires a consistent gene table across samples."
          ),
          basename(files[[1]]), basename(files[[i]])
        )
      }
    }
  }
  ref
}

.study_format_version <- 1L

.new_study <- function(txi, metadata, files, source_dir) {
  structure(
    list(
      format_version = .study_format_version,
      rsemflow_version = .rsemflow_version(),
      created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      source_dir = normalizePath(source_dir, mustWork = TRUE),
      files = unname(normalizePath(files, mustWork = TRUE)),
      sample_ids = names(files),
      txi = txi,
      metadata = metadata,
      annotation = NULL
    ),
    class = "rsemflow_study"
  )
}

#' Import RSEM gene-level results and metadata
#'
#' Finds every `*.genes.results` file in `rsem_dir`, takes the sample ID from
#' the file name (`S01.genes.results` -> `S01`), matches the files to the
#' `sample_id` column of `metadata`, and imports the counts with
#' `tximport(type = "rsem")`.
#'
#' Import stops when a file lacks the standard RSEM columns, when sample IDs
#' do not match the metadata one-to-one, or when the files do not share an
#' identical gene table (same IDs in the same order, without duplicates).
#'
#' Metadata columns whose values all parse as numbers become numeric; the rest
#' become factors. Use `categorical`/`numeric` to override.
#'
#' @param rsem_dir Directory containing `*.genes.results` files.
#' @param metadata Path to TSV or CSV metadata with a `sample_id` column.
#' @param categorical Optional metadata columns to force to categorical.
#' @param numeric Optional metadata columns to force to numeric.
#' @return An object of class `rsemflow_study`.
#' @examples
#' ex <- system.file("extdata", "example", package = "rsemflow")
#' study <- read_rsem_study(file.path(ex, "data"), file.path(ex, "metadata.tsv"))
#' study
#' @export
read_rsem_study <- function(rsem_dir, metadata, categorical = character(), numeric = character()) {
  files <- .find_rsem_files(rsem_dir)
  invisible(lapply(files, .validate_rsem_header))

  md <- .read_metadata(metadata)
  file_ids <- names(files)

  missing_md <- setdiff(file_ids, md$sample_id)
  missing_files <- setdiff(md$sample_id, file_ids)
  if (length(missing_md) || length(missing_files)) {
    msg <- c()
    if (length(missing_md)) {
      msg <- c(msg, sprintf("RSEM samples missing from metadata: %s", paste(missing_md, collapse = ", ")))
    }
    if (length(missing_files)) {
      msg <- c(msg, sprintf("Metadata samples missing RSEM files: %s", paste(missing_files, collapse = ", ")))
    }
    .stopf("%s", paste(msg, collapse = "\n"))
  }

  md <- md[match(file_ids, md$sample_id), , drop = FALSE]
  md <- .coerce_metadata(md, categorical = categorical, numeric = numeric)
  rownames(md) <- md$sample_id

  .validate_rsem_gene_ids(files)

  txi <- .quietly(tximport::tximport(
    files = files,
    type = "rsem",
    txIn = FALSE,
    txOut = FALSE
  ))

  if (!identical(colnames(txi$counts), file_ids)) {
    colnames(txi$counts) <- file_ids
    colnames(txi$abundance) <- file_ids
    colnames(txi$length) <- file_ids
  }

  .new_study(txi, md, files, rsem_dir)
}

#' Read and write RSEMflow studies
#'
#' `write_study()` saves the study as `study.rds` plus plain TSV copies of the
#' metadata, the expected-count, TPM, and effective-length matrices, a
#' sample-to-file table, and `study-info.tsv`. `read_study()` restores the
#' `study.rds`.
#'
#' @param path Study directory or a `study.rds` path.
#' @return `read_study()` returns an object of class `rsemflow_study`;
#'   `write_study()` invisibly returns `path`.
#' @examples
#' ex <- system.file("extdata", "example", package = "rsemflow")
#' study <- read_rsem_study(file.path(ex, "data"), file.path(ex, "metadata.tsv"))
#' out <- file.path(tempdir(), "study")
#' write_study(study, out)
#' list.files(out)
#' identical(read_study(out)$sample_ids, study$sample_ids)
#' @export
read_study <- function(path) {
  rds <- if (dir.exists(path)) file.path(path, "study.rds") else path
  if (!file.exists(rds)) {
    .stopf("Study file does not exist: %s", rds)
  }
  x <- readRDS(rds)
  if (!inherits(x, "rsemflow_study")) {
    .stopf("Input is not an RSEMflow study: %s", rds)
  }
  if (!identical(as.integer(x$format_version), .study_format_version)) {
    .stopf(
      "Study format version %s is not supported by rsemflow %s (expected %d). Re-run 'rsemflow data import'.",
      x$format_version, .rsemflow_version(), .study_format_version
    )
  }
  x
}

#' @param study A `rsemflow_study`.
#' @rdname read_study
#' @export
write_study <- function(study, path) {
  if (!inherits(study, "rsemflow_study")) {
    .stopf("study must be an rsemflow_study object.")
  }
  .ensure_dir(path)
  saveRDS(study, file.path(path, "study.rds"), version = 3)
  .write_tsv(study$metadata, file.path(path, "metadata.tsv"))
  .write_tsv(
    data.frame(
      sample_id = study$sample_ids,
      rsem_file = basename(study$files),
      stringsAsFactors = FALSE
    ),
    file.path(path, "samples.tsv")
  )
  .write_tsv(.matrix_to_table(study$txi$counts), file.path(path, "expected_count.tsv"))
  .write_tsv(.matrix_to_table(study$txi$abundance), file.path(path, "tpm.tsv"))
  .write_tsv(.matrix_to_table(study$txi$length), file.path(path, "effective_length.tsv"))
  .write_key_value(
    list(
      rsemflow_version = study$rsemflow_version,
      format_version = study$format_version,
      created_at = study$created_at,
      samples = length(study$sample_ids),
      genes = nrow(study$txi$counts),
      zero_effective_length_values = sum(study$txi$length <= 0, na.rm = TRUE)
    ),
    file.path(path, "study-info.tsv")
  )
  invisible(path)
}

#' @export
print.rsemflow_study <- function(x, ...) {
  md_cols <- setdiff(names(x$metadata), "sample_id")
  cat("<rsemflow_study>\n")
  cat(sprintf("  samples:  %d\n", length(x$sample_ids)))
  cat(sprintf("  genes:    %d\n", nrow(x$txi$counts)))
  cat(sprintf("  metadata: %s\n", if (length(md_cols)) paste(md_cols, collapse = ", ") else "(none)"))
  cat(sprintf("  created:  %s (rsemflow %s)\n", x$created_at, x$rsemflow_version))
  invisible(x)
}
