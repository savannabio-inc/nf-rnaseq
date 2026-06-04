#!/usr/bin/env bash
# Run nf-test using the project-local .venv-test (Nextflow + nf-test; no conda).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

VENV="${ROOT}/.venv-test"
if [[ ! -x "${VENV}/bin/nf-test" || ! -x "${VENV}/bin/nextflow" ]]; then
  echo "Test venv not set up. Running setup-test-venv.sh..."
  "${ROOT}/scripts/setup-test-venv.sh"
fi

# Prefer project-local Nextflow and nf-test over anything else on PATH
export PATH="${VENV}/bin:${PATH}"
export NFT_WORKDIR="${NFT_WORKDIR:-${ROOT}/.nf-test}"

ARGS=("$@")
if [[ ${#ARGS[@]} -eq 0 ]]; then
  ARGS=(
    "subworkflows/local/utils_nfcore_rnaseq_pipeline/tests/direct_fastq_input.nf.test"
  )
fi

echo "nf-test $(nf-test version 2>/dev/null | head -1 || true)"
echo "nextflow $(nextflow -version 2>&1 | head -1)"
echo "Running: nf-test test --verbose ${ARGS[*]}"
nf-test test --verbose "${ARGS[@]}"
