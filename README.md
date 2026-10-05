# RSEMflow

**RSEMflow** is a local command-line workflow for downstream analysis of RSEM
`*.genes.results` files using ordinary TSV/CSV sample metadata.

It is designed for experimental researchers who want standard bulk RNA-seq
downstream calculations without writing analysis scripts and without generating
figures. All outputs are TSV tables.

## What it computes

- RSEM import and input validation
- Ensembl ID annotation (`ENSG...xx` -> stable Ensembl ID + gene symbol)
- expression-level sample and gene summaries
- normalized counts and variance-stabilized expression
- PCA scores, loadings, and explained variance
- sample correlation
- DESeq2 differential expression from metadata-driven formulas
- preranked GSEA using DESeq2 statistics
- GSVA or ssGSEA sample-level pathway scores
- MSigDB access through `msigdbr`

## Local installation

RSEMflow is intended to be installed locally from its source directory.

```bash
cd rsemflow
bash install_local.sh
```

The installer:

1. checks that R is available;
2. installs missing CRAN/Bioconductor dependencies;
3. runs `R CMD INSTALL .`;
4. links the `rsemflow` executable into `${RSEMFLOW_BIN:-$HOME/.local/bin}`.

If `$HOME/.local/bin` is not in your `PATH`:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Then:

```bash
rsemflow --help
```

## Input

```text
project/
├── data/
│   ├── S01.genes.results
│   ├── S02.genes.results
│   ├── S03.genes.results
│   └── S04.genes.results
└── metadata.tsv
```

Example metadata:

```text
sample_id	condition	batch	age
S01	Control	A	51
S02	Control	B	43
S03	Drug	A	48
S04	Drug	B	55
```

No special metadata header or QIIME2-specific syntax is required.

RSEMflow infers numeric columns as numeric and other columns as categorical.
If a categorical variable is encoded with numbers (for example `batch = 1,2,3`
or `time = 0,6,24` when those are discrete time points), explicitly mark it:

```bash
--p-categorical batch --p-categorical time
```

## CLI grammar

RSEMflow uses a consistent action grammar:

```text
rsemflow <module> <action> [options]

--i-*   input data
--m-*   metadata
--p-*   analysis parameters
--o-*   outputs
```

Main actions:

```text
rsemflow data import
rsemflow data summarize

rsemflow annotation ensembl
rsemflow annotation gtf

rsemflow expression normalize
rsemflow expression pca
rsemflow expression correlation

rsemflow differential deseq2

rsemflow enrichment gsea
rsemflow enrichment gsva
```

## 1. Import RSEM data

```bash
rsemflow data import \
  --i-rsem-dir data/ \
  --m-metadata metadata.tsv \
  --o-study study/
```

This validates sample matching and RSEM format, imports RSEM with `tximport`,
and writes:

```text
study/
├── study.rds
├── metadata.tsv
├── samples.tsv
├── expected_count.tsv
├── tpm.tsv
└── effective_length.tsv
```

## 2. Expression summaries

```bash
rsemflow data summarize \
  --i-study study/ \
  --o-summary results/summary/
```

## 3. Annotation

Human Ensembl annotation:

```bash
rsemflow annotation ensembl \
  --i-study study/ \
  --p-species human \
  --o-annotation results/annotation.tsv
```

Using the same GTF that was used to build the RSEM reference is preferable
when available:

```bash
rsemflow annotation gtf \
  --i-study study/ \
  --i-gtf Homo_sapiens.GRCh38.gtf.gz \
  --o-annotation results/annotation.tsv
```

RSEMflow never replaces the original gene ID. It keeps:

- `gene_id`: original RSEM gene ID
- `ensembl_gene_id`: version suffix removed
- `gene_symbol`: annotation
- additional annotation fields when available

## 4. PCA and correlation

```bash
rsemflow expression pca \
  --i-study study/ \
  --o-pca results/pca/
```

Default PCA uses the 500 most variable genes after DESeq2 variance-stabilizing
transformation.

```bash
rsemflow expression correlation \
  --i-study study/ \
  --o-correlation results/correlation.tsv
```

Default correlation is Pearson correlation on variance-stabilized expression.

## 5. Differential expression

Simple two-group experiment:

```bash
rsemflow differential deseq2 \
  --i-study study/ \
  --p-formula '~ condition' \
  --p-term condition \
  --p-reference condition::Control \
  --o-differential results/differential/
```

Batch-adjusted analysis:

```bash
rsemflow differential deseq2 \
  --i-study study/ \
  --p-formula '~ batch + condition' \
  --p-term condition \
  --p-reference condition::Control \
  --o-differential results/differential/
```

Paired experiment:

```bash
rsemflow differential deseq2 \
  --i-study study/ \
  --p-formula '~ patient + condition' \
  --p-term condition \
  --p-reference condition::Pre \
  --o-differential results/differential/
```

Continuous variable:

```bash
rsemflow differential deseq2 \
  --i-study study/ \
  --p-formula '~ sex + age + response_score' \
  --p-term response_score \
  --o-differential results/differential/
```

Three groups, with an additional explicit comparison:

```bash
rsemflow differential deseq2 \
  --i-study study/ \
  --p-formula '~ treatment' \
  --p-term treatment \
  --p-reference treatment::Vehicle \
  --p-contrast treatment::DrugB::DrugA \
  --o-differential results/differential/
```

Interaction:

```bash
rsemflow differential deseq2 \
  --i-study study/ \
  --p-formula '~ genotype * treatment' \
  --p-term genotype:treatment \
  --p-reference genotype::WT \
  --p-reference treatment::Vehicle \
  --o-differential results/differential/
```

LRT / omnibus test:

```bash
rsemflow differential deseq2 \
  --i-study study/ \
  --p-formula '~ patient + time' \
  --p-test lrt \
  --p-reduced '~ patient' \
  --p-term time \
  --o-differential results/differential/
```

## 6. GSEA

```bash
rsemflow enrichment gsea \
  --i-differential results/differential/condition_Drug_vs_Control.tsv \
  --p-species human \
  --p-collection hallmark \
  --o-enrichment results/gsea.tsv
```

Built-in collection aliases:

- `hallmark`
- `reactome`
- `go-bp`
- `go-mf`
- `go-cc`

A raw MSigDB collection can also be used, for example:

```bash
--p-collection C2 --p-subcollection CP:REACTOME
```

Custom GMT files are accepted with `--i-genesets`. `msigdbr` may download/cache the selected MSigDB database on first use; custom GMT analysis can be run without that download.

## 7. GSVA / ssGSEA

```bash
rsemflow enrichment gsva \
  --i-study study/ \
  --p-species human \
  --p-collection hallmark \
  --p-method gsva \
  --o-scores results/gsva.tsv
```

Or:

```bash
--p-method ssgsea
```

The default expression input is `log2(TPM + 1)`.

## Help

```bash
rsemflow --help
rsemflow data --help
rsemflow data import --help
rsemflow differential deseq2 --help
```

## R API

The command-line interface is a thin wrapper around ordinary R functions:

```r
library(rsemflow)

study <- read_rsem_study("data", "metadata.tsv")
write_study(study, "study")

compute_pca(study)
compute_correlation(study)

de <- differential_deseq2(
    study,
    formula = ~ batch + condition,
    term = "condition",
    reference = "condition::Control"
)
```

## Notes

- RSEMflow performs expression-level QC, not alignment/read-level QC.
- It does not automatically remove samples.
- It does not automatically choose the experimental design.
- The user specifies the statistical model through metadata and a formula.
- No plots, PDF files, HTML reports, or interactive interfaces are produced.
