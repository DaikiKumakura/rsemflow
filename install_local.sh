#!/usr/bin/env bash
set -euo pipefail

if ! command -v R >/dev/null 2>&1 || ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: R and Rscript must be installed before installing RSEMflow." >&2
  exit 1
fi

PKG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${RSEMFLOW_BIN:-$HOME/.local/bin}"

echo "[1/4] Installing missing R/Bioconductor dependencies..."
Rscript "$PKG_DIR/scripts/install_dependencies.R"

echo "[2/4] Installing rsemflow from local source..."
R CMD INSTALL "$PKG_DIR"

echo "[3/4] Installing CLI launcher..."
mkdir -p "$BIN_DIR"
EXEC_PATH="$(Rscript --vanilla -e 'cat(system.file("exec", "rsemflow", package = "rsemflow"))')"
if [[ -z "$EXEC_PATH" || ! -f "$EXEC_PATH" ]]; then
  echo "ERROR: Installed rsemflow executable was not found." >&2
  exit 1
fi
chmod +x "$EXEC_PATH"
ln -sf "$EXEC_PATH" "$BIN_DIR/rsemflow"

echo "[4/4] Verifying installation..."
Rscript --vanilla -e 'library(rsemflow); cat("rsemflow ", as.character(packageVersion("rsemflow")), "\n", sep="")'

echo
echo "Installed CLI: $BIN_DIR/rsemflow"
case ":$PATH:" in
  *":$BIN_DIR:"*) ;;
  *)
    echo "Add this directory to PATH if needed:"
    echo "  export PATH=\"$BIN_DIR:\$PATH\""
    ;;
esac
echo
echo "Try:"
echo "  rsemflow --help"
