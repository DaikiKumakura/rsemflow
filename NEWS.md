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
