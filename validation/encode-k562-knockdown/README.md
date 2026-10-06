# Real-data validation: ENCODE K562 shRNA knockdowns

This check runs rsemflow on real RSEM output where part of the answer is
known: a knockdown should reduce its own target gene.

## Data

Six RSEM gene quantification files from ENCODE (human K562 cells, GRCh38,
GENCODE V29), in the same column layout as RSEM `*.genes.results`:

| Group | ENCODE experiment | Files (rep 1, rep 2) |
| --- | --- | --- |
| Control (non-targeting shRNA) | ENCSR129RWD | ENCFF587DSG, ENCFF774RTC |
| PTBP1 knockdown | ENCSR527IVX | ENCFF075OME, ENCFF896CIS |
| SRSF1 knockdown | ENCSR066VOO | ENCFF711WSO, ENCFF322CJL |

Both knockdown experiments list ENCSR129RWD as their control. Accessions and
MD5 sums are in [`samples.tsv`](samples.tsv). The files (about 65 MB) are not
stored in this repository; `run_validation.sh` downloads them and checks the
MD5 sums.

Each file has 59,526 rows: 58,735 versioned Ensembl gene IDs, 97 ERCC
spike-ins, 649 numeric IDs, and 45 `_PAR_Y` IDs. 35,592 effective-length
values are 0.

## Run

```bash
bash validation/encode-k562-knockdown/run_validation.sh encode-validation
```

The script runs import, summaries, PCA, correlation, a three-group DESeq2
analysis (`~ condition`, reference Control, plus an explicit SRSF1_KD vs
PTBP1_KD contrast), OrgDb annotation, hallmark GSEA for each knockdown, and
hallmark GSVA, then runs [`check_results.R`](check_results.R).

## Result (October 6, 2026)

rsemflow 0.1.1, R 4.5.3 on Windows. All 13 checks passed.

| Check | Result |
| --- | --- |
| PTBP1 in PTBP1 knockdown vs control | log2FC −1.85, padj 8.2e-112, 6th smallest p-value of 59,526 |
| SRSF1 in SRSF1 knockdown vs control | log2FC −2.55, padj 1.4e-207, 2nd smallest p-value |
| Each target in the other knockdown | SRSF1 +0.28 in PTBP1 knockdown; PTBP1 −0.20 in SRSF1 knockdown |
| PTBP2 in PTBP1 knockdown (PTBP1 represses PTBP2) | log2FC +2.54, padj 1.2e-94 |
| Explicit contrast SRSF1_KD vs PTBP1_KD | SRSF1 −2.84, PTBP1 +1.65 |
| PCA (PC1 52%, PC2 39% of variance) | every sample's nearest neighbour is its replicate |
| Annotation | 36,385 mapped; 23,141 unmapped, including all 791 non-Ensembl rows |
| Hallmark GSEA / GSVA | 50 of 50 sets; 89.6% of hallmark genes are in the ranked list |

Each command took 3–19 seconds.

These are checks of the software against known effects, not a biological
analysis: there are two replicates per group, and the PCA and GSEA results are
not interpreted here.

## Problems this check found

Fixed in rsemflow 0.1.1:

- `annotation ensembl` reported 22,350 genes as `mapped` although the OrgDb had
  no symbol, Entrez ID, or description for them.
- On Windows, `enrichment gsea` printed a progress bar and started a worker
  process.
- `annotation ensembl` printed an empty line when loading the OrgDb package.
