# rsemflow 0.1.1

Fixes found by running rsemflow on real ENCODE RSEM output
(`validation/encode-k562-knockdown/`).

* `annotation ensembl`: genes are `mapped` only when the OrgDb returned a
  symbol, Entrez ID, or description. Previously, Ensembl IDs absent from the
  OrgDb were also reported as `mapped` (22,350 of 59,526 genes in the
  validation data).
* `enrichment gsea` runs serially without a progress bar. On Windows it
  previously printed a progress bar and started a worker process. Results are
  unchanged.
* `annotation ensembl` no longer prints an empty line when loading the OrgDb.
* Added a reproducible real-data validation on ENCODE K562 PTBP1 and SRSF1
  knockdowns.

# rsemflow 0.1.0

First public release.

* `data import`: RSEM `*.genes.results` import through `tximport`, with checks
  for required columns, one-to-one sample matching, and identical gene tables.
  `sample_id` is always read as text.
* `data summarize`: sample- and gene-level expression summaries.
* `annotation ensembl` / `annotation gtf`: gene annotation that keeps the
  original RSEM gene ID.
* `expression normalize` / `pca` / `correlation`: VST, size-factor normalized
  counts, or `log2(TPM + 1)`; PCA tables; sample correlation.
* `differential deseq2`: metadata-driven DESeq2 with explicit formula, term,
  and reference levels; multi-group, numeric, interaction, simple-effect
  (`--p-at`), and LRT analyses. Coefficients are resolved from the design
  matrix, and full-rank designs are required.
* `enrichment gsea` / `gsva`: preranked GSEA (`fgseaMultilevel`, fixed seed)
  and GSVA/ssGSEA scores with MSigDB (`msigdbr`, species-aware collection
  aliases) or custom GMT gene sets.
* Genes with RSEM `effective_length = 0` are handled when building DESeq2
  objects.
* Supports human and mouse data.
* All results are TSV tables; no figures are produced.
