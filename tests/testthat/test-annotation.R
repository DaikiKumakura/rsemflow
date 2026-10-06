test_that("GTF annotation matches versioned and stable IDs", {
  study <- example_study()
  gtf <- tempfile(fileext = ".gtf")
  gtf_line <- function(type, attrs) {
    paste("1", "src", type, "1", "100", ".", "+", ".", attrs, sep = "\t")
  }
  writeLines(c(
    "#!genome-build GRCh38",
    gtf_line("gene", 'gene_id "ENSG00000000001.2"; gene_name "GENE1"; gene_type "protein_coding";'),
    gtf_line("transcript", 'gene_id "ENSG00000000001.2"; transcript_id "ENST00000000001.1";'),
    gtf_line("gene", 'gene_id "ENSG00000000002"; gene_name "GENE2"; gene_biotype "lncRNA";')
  ), gtf)
  ann <- annotate_gtf(study, gtf)
  expect_equal(nrow(ann), 200L)
  expect_equal(ann$gene_id, rownames(study$txi$counts))
  expect_equal(ann$gene_symbol[1:2], c("GENE1", "GENE2"))
  expect_equal(ann$gene_biotype[1:2], c("protein_coding", "lncRNA"))
  expect_equal(sum(ann$annotation_status == "mapped"), 2L)
})

test_that("OrgDb annotation keeps the original gene IDs", {
  skip_if_not_installed("org.Hs.eg.db")
  study <- example_study()
  ann <- annotate_ensembl(study, species = "human")
  expect_equal(ann$gene_id, rownames(study$txi$counts))
  expect_equal(ann$ensembl_gene_id, rsemflow:::.strip_ensembl_version(ann$gene_id))
  has_info <- !is.na(ann$gene_symbol) | !is.na(ann$entrez_gene_id) | !is.na(ann$gene_description)
  expect_equal(ann$annotation_status == "mapped", has_info)
  # The synthetic IDs include ones that do not exist in Ensembl.
  expect_true(any(ann$annotation_status == "unmapped"))
  expect_equal(unique(ann$annotation_source), "org.Hs.eg.db")
})
