example_dir <- function() {
  system.file("extdata", "example", package = "rsemflow")
}

example_study <- function() {
  ex <- example_dir()
  read_rsem_study(file.path(ex, "data"), file.path(ex, "metadata.tsv"))
}

# Write a minimal RSEM gene-level dataset: one *.genes.results file per column
# of `counts` and a metadata TSV. Returns the directory.
write_rsem_dataset <- function(counts, metadata, effective_length = NULL, dir = tempfile("rsem")) {
  dir.create(file.path(dir, "data"), recursive = TRUE)
  gene_ids <- rownames(counts)
  if (is.null(effective_length)) effective_length <- rep(1000, nrow(counts))
  for (s in colnames(counts)) {
    cnt <- counts[, s]
    rate <- ifelse(effective_length > 0, cnt / pmax(effective_length, 1), 0)
    tpm <- rate / sum(rate) * 1e6
    tab <- data.frame(
      gene_id = gene_ids,
      `transcript_id(s)` = sub("G", "T", gene_ids),
      length = effective_length + 150,
      effective_length = effective_length,
      expected_count = cnt,
      TPM = tpm,
      FPKM = tpm,
      check.names = FALSE
    )
    utils::write.table(
      tab, file.path(dir, "data", paste0(s, ".genes.results")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  utils::write.table(
    metadata, file.path(dir, "metadata.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  dir
}

simulate_counts <- function(n_genes, sample_ids, log2_effect = NULL, seed = 1L) {
  set.seed(seed)
  base <- stats::rgamma(n_genes, shape = 2, rate = 0.01) + 20
  mu <- matrix(base, n_genes, length(sample_ids))
  if (!is.null(log2_effect)) mu <- mu * 2^log2_effect
  counts <- matrix(
    stats::rnbinom(length(mu), mu = mu, size = 20),
    n_genes, length(sample_ids),
    dimnames = list(sprintf("ENSG%011d.1", seq_len(n_genes)), sample_ids)
  )
  counts
}

# 2 x 2 genotype-by-treatment design with 3 replicates per cell. The first 50
# genes respond to Drug only in KO samples.
interaction_dataset <- function() {
  md <- expand.grid(rep = 1:3, treatment = c("Vehicle", "Drug"), genotype = c("WT", "KO"))
  md$sample_id <- sprintf("%s_%s_%d", md$genotype, md$treatment, md$rep)
  md <- md[, c("sample_id", "genotype", "treatment")]
  effect <- matrix(0, 300, nrow(md))
  effect[1:50, md$genotype == "KO" & md$treatment == "Drug"] <- 2
  counts <- simulate_counts(300, md$sample_id, effect, seed = 11L)
  write_rsem_dataset(counts, md)
}
