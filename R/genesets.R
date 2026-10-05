.collection_alias <- function(collection, subcollection = NULL) {
  key <- tolower(trimws(collection))
  aliases <- list(
    hallmark = c("H", ""),
    reactome = c("C2", "CP:REACTOME"),
    "go-bp" = c("C5", "GO:BP"),
    "go-mf" = c("C5", "GO:MF"),
    "go-cc" = c("C5", "GO:CC")
  )
  if (key %in% names(aliases)) {
    return(list(collection = aliases[[key]][[1]], subcollection = aliases[[key]][[2]]))
  }
  list(collection = collection, subcollection = subcollection %||% "")
}

.read_gmt <- function(path) {
  if (!file.exists(path)) .stopf("GMT file does not exist: %s", path)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  sets <- lapply(strsplit(lines, "\t", fixed = TRUE), function(x) {
    if (length(x) < 3L) return(NULL)
    unique(x[-c(1, 2)])
  })
  nms <- vapply(strsplit(lines, "\t", fixed = TRUE), function(x) x[[1]], character(1))
  keep <- !vapply(sets, is.null, logical(1))
  sets <- sets[keep]
  names(sets) <- nms[keep]
  if (!length(sets)) .stopf("No gene sets were read from GMT: %s", path)
  sets
}

.get_msigdbr_genesets <- function(species, collection = "hallmark", subcollection = NULL, ids = NULL) {
  species <- .resolve_species(species, ids)
  resolved <- .collection_alias(collection, subcollection)
  scientific <- .species_scientific(species)
  db_species <- .species_db(species)

  args <- list(
    db_species = db_species,
    species = scientific,
    collection = resolved$collection
  )
  if (nzchar(resolved$subcollection)) {
    args$subcollection <- resolved$subcollection
  }

  tbl <- do.call(msigdbr::msigdbr, args)
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  if (!nrow(tbl)) {
    .stopf(
      "No MSigDB gene sets were returned for collection '%s' and subcollection '%s'.",
      resolved$collection, resolved$subcollection
    )
  }
  tbl <- tbl[!is.na(tbl$ensembl_gene) & nzchar(tbl$ensembl_gene), , drop = FALSE]
  sets <- split(tbl$ensembl_gene, tbl$gs_name)
  sets <- lapply(sets, unique)

  info <- unique(tbl[, intersect(
    c("gs_name", "gs_collection", "gs_subcollection", "gs_collection_name", "gs_description", "db_version", "db_target_species"),
    names(tbl)
  ), drop = FALSE])
  list(
    sets = sets,
    info = info,
    collection = resolved$collection,
    subcollection = resolved$subcollection,
    db_version = unique(tbl$db_version)
  )
}

.resolve_genesets <- function(
    genesets = NULL,
    species = "auto",
    collection = "hallmark",
    subcollection = NULL,
    ids = NULL) {
  if (!is.null(genesets)) {
    sets <- if (is.character(genesets) && length(genesets) == 1L) .read_gmt(genesets) else genesets
    if (!is.list(sets) || is.null(names(sets))) {
      .stopf("Custom gene sets must be a named list or a GMT path.")
    }
    return(list(
      sets = sets,
      info = data.frame(gs_name = names(sets), stringsAsFactors = FALSE),
      collection = "custom",
      subcollection = "",
      db_version = NA_character_
    ))
  }
  .get_msigdbr_genesets(species, collection, subcollection, ids = ids)
}
