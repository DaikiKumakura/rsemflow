# Checks for the ENCODE K562 knockdown validation. Run from run_validation.sh,
# or as `Rscript check_results.R <work-directory>`.
#
# Expected results, from known biology and from how the data are organized:
# * each shRNA target is strongly reduced in its own knockdown and among the
#   most significant genes there, and changes little in the other knockdown;
# * PTBP2 rises when PTBP1 is knocked down (PTBP1 represses PTBP2);
# * PCA separates the three groups, with replicates closest to each other;
# * non-Ensembl rows (spike-ins, numeric IDs, _PAR_Y) are unmapped, and a
#   mapped gene always has some annotation;
# * most hallmark genes are present in the ranked list.

args <- commandArgs(trailingOnly = TRUE)
work <- if (length(args)) args[[1]] else "."
path <- function(...) file.path(work, ...)

ids <- list(PTBP1 = "ENSG00000011304", SRSF1 = "ENSG00000136450", PTBP2 = "ENSG00000117569")
de_table <- function(name) {
  utils::read.delim(path("results", "differential", paste0(name, ".tsv")))
}
gene_row <- function(de, gene) de[de$ensembl_gene_id == ids[[gene]], ]

checks <- list()
check <- function(label, ok, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(check = label, pass = isTRUE(ok), detail = detail)
}

info <- utils::read.delim(path("study", "study-info.tsv"))
n_genes <- as.integer(info$value[info$key == "genes"])
check("all genes imported", n_genes == 59526L, sprintf("%d genes", n_genes))

ptbp1_kd <- de_table("condition_PTBP1_KD_vs_Control")
srsf1_kd <- de_table("condition_SRSF1_KD_vs_Control")
kd_vs_kd <- de_table("condition_SRSF1_KD_vs_PTBP1_KD")

target_check <- function(de, gene, label) {
  r <- gene_row(de, gene)
  rank_p <- rank(de$pvalue, na.last = "keep")[de$ensembl_gene_id == ids[[gene]]]
  check(
    label,
    r$log2FoldChange < -1 && r$padj < 1e-10 && rank_p <= 20,
    sprintf("log2FC %.2f, padj %.1e, p-value rank %d", r$log2FoldChange, r$padj, as.integer(rank_p))
  )
}
target_check(ptbp1_kd, "PTBP1", "PTBP1 down in PTBP1 knockdown")
target_check(srsf1_kd, "SRSF1", "SRSF1 down in SRSF1 knockdown")

off <- function(de, gene, label) {
  r <- gene_row(de, gene)
  check(label, abs(r$log2FoldChange) < 0.5, sprintf("log2FC %.2f", r$log2FoldChange))
}
off(ptbp1_kd, "SRSF1", "SRSF1 nearly unchanged in PTBP1 knockdown")
off(srsf1_kd, "PTBP1", "PTBP1 nearly unchanged in SRSF1 knockdown")

r <- gene_row(ptbp1_kd, "PTBP2")
check("PTBP2 up in PTBP1 knockdown", r$log2FoldChange > 1 && r$padj < 1e-10,
      sprintf("log2FC %.2f, padj %.1e", r$log2FoldChange, r$padj))

r1 <- gene_row(kd_vs_kd, "SRSF1"); r2 <- gene_row(kd_vs_kd, "PTBP1")
check("explicit contrast SRSF1_KD vs PTBP1_KD has the expected signs",
      r1$log2FoldChange < -1 && r2$log2FoldChange > 1,
      sprintf("SRSF1 %.2f, PTBP1 %.2f", r1$log2FoldChange, r2$log2FoldChange))

scores <- utils::read.delim(path("results", "pca", "scores.tsv"))
xy <- as.matrix(scores[, c("PC1", "PC2")]); rownames(xy) <- scores$sample_id
d <- as.matrix(stats::dist(xy)); diag(d) <- Inf
nearest <- colnames(d)[apply(d, 1, which.min)]
group <- sub("_[0-9]+$", "", rownames(d))
check("each sample's nearest PCA neighbour is its replicate",
      all(group == sub("_[0-9]+$", "", nearest)),
      paste(rownames(d), nearest, sep = "->", collapse = ", "))

ann <- utils::read.delim(path("results", "annotation.tsv"))
non_ensembl <- !grepl("^ENSG[0-9]+$", ann$ensembl_gene_id)
has_info <- !is.na(ann$gene_symbol) | !is.na(ann$entrez_gene_id) | !is.na(ann$gene_description)
check("non-Ensembl rows are unmapped", all(ann$annotation_status[non_ensembl] == "unmapped"),
      sprintf("%d non-Ensembl rows", sum(non_ensembl)))
check("mapped means some annotation was found", all((ann$annotation_status == "mapped") == has_info),
      sprintf("%d mapped, %d unmapped", sum(ann$annotation_status == "mapped"), sum(ann$annotation_status == "unmapped")))
sym <- ann$gene_symbol[match(paste0(ids$PTBP1), ann$ensembl_gene_id)]
check("PTBP1 symbol", identical(sym, "PTBP1"), sym)

gsea <- utils::read.delim(path("results", "gsea_PTBP1_hallmark.tsv"))
check("GSEA returns all 50 hallmark sets", nrow(gsea) == 50L, sprintf("%d sets", nrow(gsea)))
gsva <- utils::read.delim(path("results", "gsva_hallmark.tsv"), check.names = FALSE)
check("GSVA scores complete", nrow(gsva) == 50L && !anyNA(gsva[, -(1:4)]),
      sprintf("%d x %d", nrow(gsva), ncol(gsva) - 4L))

out <- do.call(rbind, checks)
print(out, right = FALSE, row.names = FALSE)
if (!all(out$pass)) {
  stop(sum(!out$pass), " check(s) failed.", call. = FALSE)
}
cat("\nAll", nrow(out), "checks passed.\n")
