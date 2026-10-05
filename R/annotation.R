
.orgdb_package_for_species <- function(species) {
  switch(
    .resolve_species(species),
    human = "org.Hs.eg.db",
    mouse = "org.Mm.eg.db"
  )
}

.orgdb_for_species <- function(species) {
  species <- .resolve_species(species)
  pkg <- .orgdb_package_for_species(species)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    .stopf(
      "Annotation package '%s' is not installed. Install it with BiocManager::install('%s').",
      pkg, pkg
    )
  }
  getExportedValue(pkg, pkg)
}

#' Annotate Ensembl gene IDs using a local Bioconductor OrgDb
#'
#' Adds the stable Ensembl ID (version suffix removed), gene symbol, Entrez ID,
#' and description from `org.Hs.eg.db` or `org.Mm.eg.db`. The
#' original RSEM `gene_id` is always kept. When an Ensembl ID maps to several
#' Entrez records, the record with the most filled fields is kept.
#'
#' @param study A `rsemflow_study` or study path.
#' @param species `auto`, `human`, or `mouse`.
#' @return A gene annotation data frame.
#' @export
annotate_ensembl <- function(study, species = "auto") {
  if (is.character(study)) study <- read_study(study)
  gene_id <- rownames(study$txi$counts)
  stable <- .strip_ensembl_version(gene_id)
  species <- .resolve_species(species, stable)
  db <- .orgdb_for_species(species)

  keys <- unique(stable[.looks_like_ensembl(stable)])
  if (!length(keys)) {
    return(data.frame(
      gene_id = gene_id,
      ensembl_gene_id = stable,
      gene_symbol = NA_character_,
      entrez_gene_id = NA_character_,
      gene_description = NA_character_,
      annotation_status = "unmapped",
      species = species,
      annotation_source = .orgdb_package_for_species(species),
      annotation_version = as.character(utils::packageVersion(.orgdb_package_for_species(species))),
      stringsAsFactors = FALSE
    ))
  }

  ann <- .quietly(AnnotationDbi::select(
    db,
    keys = keys,
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL", "ENTREZID", "GENENAME")
  ))
  ann <- as.data.frame(ann, stringsAsFactors = FALSE)

  if (nrow(ann)) {
    ann$.score <- rowSums(!is.na(ann[, intersect(c("SYMBOL", "ENTREZID", "GENENAME"), names(ann)), drop = FALSE]))
    ann <- ann[order(ann$ENSEMBL, -ann$.score), , drop = FALSE]
    ann <- ann[!duplicated(ann$ENSEMBL), , drop = FALSE]
    ann$.score <- NULL
  }

  idx <- match(stable, ann$ENSEMBL)
  symbol <- ann$SYMBOL[idx]
  entrez <- ann$ENTREZID[idx]
  description <- ann$GENENAME[idx]
  status <- ifelse(is.na(idx), "unmapped", "mapped")

  data.frame(
    gene_id = gene_id,
    ensembl_gene_id = stable,
    gene_symbol = symbol,
    entrez_gene_id = entrez,
    gene_description = description,
    annotation_status = status,
    species = species,
    annotation_source = .orgdb_package_for_species(species),
    annotation_version = as.character(utils::packageVersion(.orgdb_package_for_species(species))),
    stringsAsFactors = FALSE
  )
}

.gtf_attr <- function(x, key) {
  pattern <- paste0("(?:^|;[[:space:]]*)", key, "[[:space:]]+\"([^\"]+)\"")
  out <- rep(NA_character_, length(x))
  m <- regexec(pattern, x, perl = TRUE)
  mm <- regmatches(x, m)
  hit <- lengths(mm) >= 2L
  out[hit] <- vapply(mm[hit], `[[`, character(1), 2L)
  out
}

.read_gtf_genes <- function(gtf) {
  if (!file.exists(gtf)) .stopf("GTF file does not exist: %s", gtf)
  con <- if (grepl("\\.gz$", gtf, ignore.case = TRUE)) gzfile(gtf, "rt") else file(gtf, "rt")
  on.exit(close(con), add = TRUE)

  attrs <- character()
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[!startsWith(lines, "#")]
    if (!length(lines)) next
    fields <- strsplit(lines, "\t", fixed = TRUE)
    is_gene <- vapply(fields, function(z) length(z) >= 9L && identical(z[[3]], "gene"), logical(1))
    if (any(is_gene)) {
      attrs <- c(attrs, vapply(fields[is_gene], `[[`, character(1), 9L))
    }
  }
  if (!length(attrs)) .stopf("No gene records were found in GTF: %s", gtf)

  gene_id <- .gtf_attr(attrs, "gene_id")
  gene_name <- .gtf_attr(attrs, "gene_name")
  gene_type <- .gtf_attr(attrs, "gene_type")
  missing_type <- is.na(gene_type)
  if (any(missing_type)) {
    gene_type[missing_type] <- .gtf_attr(attrs[missing_type], "gene_biotype")
  }
  description <- .gtf_attr(attrs, "description")

  out <- data.frame(
    gtf_gene_id = gene_id,
    ensembl_gene_id = .strip_ensembl_version(gene_id),
    gene_symbol = gene_name,
    gene_biotype = gene_type,
    gene_description = description,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$gtf_gene_id) & !duplicated(out$gtf_gene_id), , drop = FALSE]
  out
}

#' Annotate genes from a GTF file
#'
#' Reads the `gene` records of the GTF used to build the RSEM reference and
#' returns `gene_name`, `gene_type`/`gene_biotype`, and `description`. IDs are
#' matched on the full versioned ID first and on the stable Ensembl ID second.
#'
#' @param study A `rsemflow_study` or study path.
#' @param gtf Path to a GTF or GTF.GZ file.
#' @return A gene annotation data frame.
#' @export
annotate_gtf <- function(study, gtf) {
  if (is.character(study)) study <- read_study(study)
  gene_id <- rownames(study$txi$counts)
  stable <- .strip_ensembl_version(gene_id)
  ann <- .read_gtf_genes(gtf)

  idx_exact <- match(gene_id, ann$gtf_gene_id)
  idx_stable <- match(stable, ann$ensembl_gene_id)
  idx <- idx_exact
  fallback <- is.na(idx)
  idx[fallback] <- idx_stable[fallback]

  data.frame(
    gene_id = gene_id,
    ensembl_gene_id = stable,
    gene_symbol = ann$gene_symbol[idx],
    gene_biotype = ann$gene_biotype[idx],
    gene_description = ann$gene_description[idx],
    annotation_status = ifelse(is.na(idx), "unmapped", "mapped"),
    annotation_source = basename(gtf),
    annotation_version = NA_character_,
    stringsAsFactors = FALSE
  )
}
