#!/usr/bin/env bash
# Run the basic rsemflow workflow on the bundled synthetic example.
# Needs no annotation database or MSigDB download.
#
#   bash run_example.sh [output-directory]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${1:-rsemflow-example-output}"

if [[ -e "$WORK" ]]; then
  echo "ERROR: '$WORK' already exists; choose a new output directory." >&2
  exit 1
fi
mkdir -p "$WORK"

rsemflow data import \
  --i-rsem-dir "$HERE/data" \
  --m-metadata "$HERE/metadata.tsv" \
  --o-study "$WORK/study"

rsemflow data summarize \
  --i-study "$WORK/study" \
  --o-summary "$WORK/summary"

rsemflow expression pca \
  --i-study "$WORK/study" \
  --o-pca "$WORK/pca"

rsemflow expression correlation \
  --i-study "$WORK/study" \
  --o-correlation "$WORK/correlation.tsv"

rsemflow differential deseq2 \
  --i-study "$WORK/study" \
  --p-formula '~ batch + condition' \
  --p-term condition \
  --p-reference condition::Control \
  --o-differential "$WORK/differential"

echo "Example analysis complete: $WORK"
