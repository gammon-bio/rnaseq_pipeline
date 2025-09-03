#!/usr/bin/env bash
set -euo pipefail

# get_refs.sh — Download Ensembl GTF and cDNA FASTA into refs/
# Created: 2025-09-03
# Usage:
#   bash scripts/get_refs.sh \
#     [--species human|mouse] [--build GRCh38|GRCm39] [--release <ensembl_release>] \
#     [--gtf_url <url>] [--fasta_url <url>]
# Defaults: --species human, --build GRCh38, --release current (Ensembl), cdna FASTA

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REFS_DIR="${ROOT_DIR}/refs"
GTF_DIR="${REFS_DIR}/gtf"
FA_DIR="${REFS_DIR}/fa"

species="human"
build="GRCh38"
release="current"   # or integer like 110
gtf_url=""
fasta_url=""

usage() {
  echo "Usage: $0 [--species human|mouse] [--build GRCh38|GRCm39] [--release <n>|current] [--gtf_url URL] [--fasta_url URL]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --species) species="$2"; shift 2 ;;
    --build) build="$2"; shift 2 ;;
    --release) release="$2"; shift 2 ;;
    --gtf_url) gtf_url="$2"; shift 2 ;;
    --fasta_url) fasta_url="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p "${GTF_DIR}" "${FA_DIR}"

# Map species/build to Ensembl naming
lower_species=""
cap_species=""
case "${species}" in
  human)
    lower_species="homo_sapiens"
    cap_species="Homo_sapiens"
    [[ "${build}" == "GRCh38" ]] || { echo "ERROR: human build must be GRCh38" >&2; exit 1; }
    ;;
  mouse)
    lower_species="mus_musculus"
    cap_species="Mus_musculus"
    [[ "${build}" == "GRCm39" ]] || { echo "ERROR: mouse build must be GRCm39" >&2; exit 1; }
    ;;
  *)
    echo "ERROR: --species must be human or mouse" >&2; exit 1 ;;
esac

base="https://ftp.ensembl.org/pub"
rel_path="${release}"
if [[ "${release}" != "current" ]]; then
  rel_path="release-${release}"
fi

# Construct URLs if not provided explicitly
if [[ -z "${gtf_url}" ]]; then
  if [[ "${release}" == "current" ]]; then
    # file names still contain release; we won't guess — use wildcard after download rename
    gtf_url="${base}/${rel_path}/gtf/${lower_species}/${cap_species}.${build}.*.gtf.gz"
  else
    gtf_url="${base}/${rel_path}/gtf/${lower_species}/${cap_species}.${build}.${release}.gtf.gz"
  fi
fi

if [[ -z "${fasta_url}" ]]; then
  if [[ "${release}" == "current" ]]; then
    fasta_url="${base}/${rel_path}/fasta/${lower_species}/cdna/${cap_species}.${build}.cdna.all.fa.gz"
  else
    fasta_url="${base}/${rel_path}/fasta/${lower_species}/cdna/${cap_species}.${build}.cdna.all.fa.gz"
  fi
fi

echo "[get_refs] Using URLs:" >&2
echo "  GTF:   ${gtf_url}" >&2
echo "  FASTA: ${fasta_url}" >&2

download() {
  local url="$1"; local out_dir="$2"
  local fname
  fname=$(basename "$url")
  local out_path="${out_dir}/${fname}"
  echo "[get_refs] Downloading ${url}" >&2
  curl -fL -C - -o "${out_path}" "${url}"
  if [[ ! -s "${out_path}" ]]; then
    echo "ERROR: Downloaded file empty: ${out_path}" >&2
    exit 1
  fi
  echo "[get_refs] sha256(${fname}): $(shasum -a 256 "${out_path}" | awk '{print $1}')" >&2
  echo "${out_path}"
}

# Handle potential wildcard in current GTF case by attempting download and glob resolution
resolve_and_download_gtf() {
  local url="$1"
  local out_dir="$2"
  if [[ "$url" == *"*"* ]]; then
    # Try known pattern with cdn path; we cannot list directory without FTP, so try without release in name fails.
    # Fallback: attempt the canonical path; if it fails, instruct user to pass --gtf_url explicitly.
    echo "[get_refs] 'current' GTF filename contains release; if this fails, pass --gtf_url explicitly." >&2
  fi
  download "$url" "$out_dir"
}

gtf_gz="$(resolve_and_download_gtf "${gtf_url}" "${GTF_DIR}")"
fa_gz="$(download "${fasta_url}" "${FA_DIR}")"

# Gunzip to plain .gtf/.fa (keep .gz files)
gtf_plain="${GTF_DIR}/$(basename "${gtf_gz}" .gz)"
fa_plain="${FA_DIR}/$(basename "${fa_gz}" .gz)"

echo "[get_refs] Decompressing archives" >&2
gzip -dc "${gtf_gz}" > "${gtf_plain}"
gzip -dc "${fa_gz}" > "${fa_plain}"

if [[ ! -s "${gtf_plain}" || ! -s "${fa_plain}" ]]; then
  echo "ERROR: Decompressed files are empty" >&2
  exit 1
fi

# Document URLs and checksums
readme="${REFS_DIR}/README.md"
{
  echo "# References"
  echo ""
  echo "- GTF URL:   ${gtf_url}"
  echo "- FASTA URL: ${fasta_url}"
  echo ""
  echo "## SHA256"
  echo "- $(basename "${gtf_gz}"): $(shasum -a 256 "${gtf_gz}" | awk '{print $1}')"
  echo "- $(basename "${fa_gz}"):  $(shasum -a 256 "${fa_gz}" | awk '{print $1}')"
  echo "- $(basename "${gtf_plain}"): $(shasum -a 256 "${gtf_plain}" | awk '{print $1}')"
  echo "- $(basename "${fa_plain}"):  $(shasum -a 256 "${fa_plain}" | awk '{print $1}')"
} > "${readme}"

echo "[get_refs] Done. Files in: ${REFS_DIR}" >&2

