#!/usr/bin/env bash
# Real-data validation on ENCODE K562 shRNA knockdown RNA-seq.
#
# Downloads six RSEM gene quantification files (about 65 MB) listed in
# samples.tsv, checks their MD5 sums, runs the rsemflow command-line workflow,
# and checks the results with check_results.R.
#
#   bash run_validation.sh [work-directory]
#
# Requires an installed rsemflow, curl, md5sum, and org.Hs.eg.db. Set
# RSEMFLOW to the command used to call the CLI if `rsemflow` is not on PATH,
# for example RSEMFLOW='Rscript -e rsemflow::cli_main()'.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${1:-encode-validation}"
RSEMFLOW="${RSEMFLOW:-rsemflow}"
mkdir -p "$WORK/raw" "$WORK/data"
cd "$WORK"

echo "== Download and verify"
printf 'sample_id\tcondition\treplicate\n' > metadata.tsv
tail -n +2 "$HERE/samples.tsv" | while IFS=$'\t' read -r sample accession md5 group replicate; do
  file="raw/$accession.tsv"
  if [[ ! -f "$file" ]]; then
    curl -sSL -o "$file" "https://www.encodeproject.org/files/$accession/@@download/$accession.tsv"
  fi
  got="$(md5sum "$file" | cut -d' ' -f1)"
  if [[ "$got" != "$md5" ]]; then
    echo "ERROR: MD5 mismatch for $accession" >&2
    exit 1
  fi
  cp "$file" "data/$sample.genes.results"
  printf '%s\t%s\t%s\n' "$sample" "$group" "$replicate" >> metadata.tsv
  echo "  $sample <- $accession (MD5 OK)"
done

run() {
  local start=$SECONDS
  $RSEMFLOW "$@"
  echo "  [$(( SECONDS - start )) s] $1 $2"
}

echo "== Run rsemflow"
rm -rf study results
run data import --i-rsem-dir data --m-metadata metadata.tsv --o-study study
run data summarize --i-study study --o-summary results/summary
run expression pca --i-study study --o-pca results/pca
run expression correlation --i-study study --o-correlation results/correlation.tsv
run differential deseq2 --i-study study \
  --p-formula '~ condition' --p-term condition \
  --p-reference condition::Control \
  --p-contrast condition::SRSF1_KD::PTBP1_KD \
  --o-differential results/differential
run annotation ensembl --i-study study --p-species human --o-annotation results/annotation.tsv
for kd in PTBP1 SRSF1; do
  run enrichment gsea \
    --i-differential "results/differential/condition_${kd}_KD_vs_Control.tsv" \
    --p-species human --p-collection hallmark \
    --o-enrichment "results/gsea_${kd}_hallmark.tsv"
done
run enrichment gsva --i-study study --p-species human --p-collection hallmark --o-scores results/gsva_hallmark.tsv

echo "== Check results"
Rscript "$HERE/check_results.R" .
