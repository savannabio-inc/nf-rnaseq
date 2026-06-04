#!/usr/bin/env bash
# Self-contained test environment: Python venv + nf-test + Nextflow (no conda required).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

VENV="${ROOT}/.venv-test"
BIN="${VENV}/bin"
NFT_VERSION="${NFT_VERSION:-0.9.0}"
NXF_VERSION="${NXF_VERSION:-24.04.2}"

if [[ ! -d "$VENV" ]]; then
  echo "Creating Python venv at ${VENV}..."
  python3 -m venv "$VENV"
fi

"${BIN}/pip" install -q -U pip
"${BIN}/pip" install -q -r "${ROOT}/requirements-test.txt"

install_to_bin() {
  local name="$1"
  shift
  echo "Installing ${name} into ${BIN}..."
  local tmp
  tmp="$(mktemp -d)"
  (cd "$tmp" && "$@")
  mv "${tmp}/${name}" "${BIN}/${name}"
  chmod +x "${BIN}/${name}"
  rm -rf "$tmp"
}

if [[ ! -x "${BIN}/nf-test" ]]; then
  install_to_bin nf-test bash -c "curl -fsSL https://get.nf-test.com | bash -s ${NFT_VERSION}"
fi

if [[ ! -x "${BIN}/nextflow" ]]; then
  install_to_bin nextflow bash -c "curl -fsSL https://get.nextflow.io | NXF_VER=${NXF_VERSION} bash"
fi

echo ""
echo "Test environment ready (no conda required):"
echo "  source ${VENV}/bin/activate"
echo "  ./scripts/run-nf-tests.sh"
echo ""
echo "nf-test:  $("${BIN}/nf-test" version 2>/dev/null | head -1 || true)"
echo "nextflow: $("${BIN}/nextflow" -version 2>&1 | head -1 || true)"
