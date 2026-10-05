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

.new_study <- function(txi, metadata, files, source_dir) {
  structure(
    list(
      format_version = 1L,
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
#' @param rsem_dir Directory containing `*.genes.results` files.
#' @param metadata Path to TSV or CSV metadata with a `sample_id` column.
#' @param categorical Optional metadata columns to force to categorical.
#' @param numeric Optional metadata columns to force to numeric.
#' @return An object of class `rsemflow_study`.
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

  txi <- tximport::tximport(
    files = files,
    type = "rsem",
    txIn = FALSE,
    txOut = FALSE
  )

  if (!identical(colnames(txi$counts), file_ids)) {
    colnames(txi$counts) <- file_ids
    colnames(txi$abundance) <- file_ids
    colnames(txi$length) <- file_ids
  }

  .new_study(txi, md, files, rsem_dir)
}

#' Read a serialized RSEMflow study
#'
#' @param path Study directory or a `study.rds` path.
#' @return An object of class `rsemflow_study`.
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
  x
}

#' Write an RSEMflow study to a local directory
#'
#' @param study A `rsemflow_study`.
#' @param path Output directory.
#' @return Invisibly returns `path`.
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
      genes = nrow(study$txi$counts)
    ),
    file.path(path, "study-info.tsv")
  )
  invisible(path)
}
