#!/usr/bin/env bash
set -euo pipefail

# get_refs.sh — Download Ensembl GTF and cDNA FASTA into data/references/
# Created: 2025-09-03
# Usage:
#   bash scripts/get_refs.sh \
#     [--species human|mouse] [--build GRCh38|GRCm39] [--release <ensembl_release>] \
#     [--gtf_flavor plain|chr|chr_patch_hapl_scaff|abinitio|auto] \
#     [--gtf_url <url>] [--fasta_url <url>] \
#     [--decoy|--no-decoy] [--genome_url URL]
# Defaults: --species human, --build GRCh38, --release 115 (Ensembl), cdna FASTA,
#           --gtf_flavor auto (prefers plain > chr_patch_hapl_scaff > chr > abinitio),
#           decoy-aware index on by default (--decoy)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REFS_DIR="${ROOT_DIR}/data/references"
GTF_DIR="${REFS_DIR}/gtf"
FA_DIR="${REFS_DIR}/fa"

species="human"
build="GRCh38"
release="115"   # or integer like 110
gtf_url=""
fasta_url=""
gtf_flavor="auto"   # plain|chr|chr_patch_hapl_scaff|abinitio|auto
decoy="true"
genome_url=""

usage() {
  echo "Usage: $0 [--species human|mouse] [--build GRCh38|GRCm39] [--release <n>|current] [--gtf_flavor FLAVOR] [--gtf_url URL] [--fasta_url URL] [--decoy|--no-decoy] [--genome_url URL]" >&2
  echo "Defaults: --release 115, decoy-aware index on (--decoy)" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --species) species="$2"; shift 2 ;;
    --build) build="$2"; shift 2 ;;
    --release) release="$2"; shift 2 ;;
    --gtf_url) gtf_url="$2"; shift 2 ;;
    --fasta_url) fasta_url="$2"; shift 2 ;;
    --gtf_flavor) gtf_flavor="$2"; shift 2 ;;
    --decoy) decoy="true"; shift 1 ;;
    --no-decoy) decoy="false"; shift 1 ;;
    --genome_url) genome_url="$2"; shift 2 ;;
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
  # Query CHECKSUMS to resolve the exact GTF filename (handles optional suffixes)
  checksums_url="${base}/${rel_path}/gtf/${lower_species}/CHECKSUMS"
  echo "[get_refs] Resolving GTF via CHECKSUMS: ${checksums_url}" >&2
  tmp_checksums="$(mktemp)"
  curl -fL -C - -o "${tmp_checksums}" "${checksums_url}"
  cand_list=$(awk '{print $NF}' "${tmp_checksums}" | grep -E "^${cap_species}\\.${build}\\.[0-9]+(\\.(chr|chr_patch_hapl_scaff|abinitio))?\\.gtf\\.gz$") || true
  rm -f "${tmp_checksums}"
  if [[ -z "${cand_list}" ]]; then
    echo "ERROR: Could not find any matching GTF in CHECKSUMS for ${cap_species}.${build}. Pass --gtf_url explicitly." >&2
    exit 1
  fi

  # Choose best match by flavor preference
  resolve_by_suffix() {
    local suffix="$1"  # e.g., "", ".chr", ".chr_patch_hapl_scaff", ".abinitio"
    local pat
    if [[ -z "$suffix" ]]; then
      # Match base GTF without any flavor suffix
      pat="^${cap_species}\\.${build}\\.[0-9]+\\.gtf\\.gz$"
    else
      # escape dots in suffix and handle versioned files
      local esc_suffix
      esc_suffix=$(printf '%s' "$suffix" | sed 's/[.]/\\./g')
      pat="^${cap_species}\\.${build}\\.[0-9]+${esc_suffix}\\.gtf\\.gz$"
    fi
    printf "%s\n" "$cand_list" | grep -E "$pat" | sort -V | tail -n1 || true
  }

  case "${gtf_flavor}" in
    auto)
      for suf in "" ".chr_patch_hapl_scaff" ".chr" ".abinitio"; do
        gtf_file_name=$(resolve_by_suffix "$suf")
        [[ -n "${gtf_file_name}" ]] && break
      done
      ;;
    plain) gtf_file_name=$(resolve_by_suffix "") ;;
    chr) gtf_file_name=$(resolve_by_suffix ".chr") ;;
    chr_patch_hapl_scaff) gtf_file_name=$(resolve_by_suffix ".chr_patch_hapl_scaff") ;;
    abinitio) gtf_file_name=$(resolve_by_suffix ".abinitio") ;;
    *) echo "ERROR: Invalid --gtf_flavor: ${gtf_flavor}" >&2; exit 1 ;;
  esac

  if [[ -z "${gtf_file_name}" ]]; then
    echo "ERROR: Could not resolve a GTF file for flavor '${gtf_flavor}'. Available candidates:" >&2
    printf "  %s\n" $cand_list >&2
    echo "Tip: pass --gtf_url to specify the exact file." >&2
    exit 1
  fi

  echo "[get_refs] Selected GTF: ${gtf_file_name} (flavor=${gtf_flavor})" >&2
  gtf_url="${base}/${rel_path}/gtf/${lower_species}/${gtf_file_name}"
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

gtf_gz="$(download "${gtf_url}" "${GTF_DIR}")"
fa_gz="$(download "${fasta_url}" "${FA_DIR}")"

# --- Decoy genome (for decoy-aware Salmon index) ---
genome_gz=""
decoys_txt=""
gentrome_gz=""
if [[ "${decoy}" == "true" ]]; then
  if [[ -z "${genome_url}" ]]; then
    dna_checksums_url="${base}/${rel_path}/fasta/${lower_species}/dna/CHECKSUMS"
    echo "[get_refs] Resolving decoy genome via CHECKSUMS: ${dna_checksums_url}" >&2
    tmp_dna="$(mktemp)"
    curl -fL -C - -o "${tmp_dna}" "${dna_checksums_url}"
    genome_name=""
    for suf in "dna.primary_assembly" "dna.toplevel"; do
      genome_name=$(awk '{print $NF}' "${tmp_dna}" | grep -E "^${cap_species}\\.${build}\\.${suf}\\.fa\\.gz$" | head -n1 || true)
      [[ -n "${genome_name}" ]] && break
    done
    rm -f "${tmp_dna}"
    if [[ -z "${genome_name}" ]]; then
      echo "ERROR: Could not resolve a primary_assembly/toplevel genome FASTA in CHECKSUMS. Pass --genome_url or use --no-decoy." >&2
      exit 1
    fi
    echo "[get_refs] Selected decoy genome: ${genome_name}" >&2
    genome_url="${base}/${rel_path}/fasta/${lower_species}/dna/${genome_name}"
  fi
  echo "[get_refs] GENOME (decoy): ${genome_url}" >&2
  genome_gz="$(download "${genome_url}" "${FA_DIR}")"

  # Build decoys.txt (genome sequence names) and gentrome.fa.gz (transcriptome + genome)
  decoys_txt="${FA_DIR}/decoys.txt"
  gentrome_gz="${FA_DIR}/gentrome.fa.gz"
  echo "[get_refs] Building decoys.txt and gentrome.fa.gz" >&2
  grep '^>' <(gzip -dc "${genome_gz}") | cut -d ' ' -f1 | sed 's/^>//' > "${decoys_txt}"
  cat "${fa_gz}" "${genome_gz}" > "${gentrome_gz}"
  if [[ ! -s "${decoys_txt}" || ! -s "${gentrome_gz}" ]]; then
    echo "ERROR: Failed to build decoys.txt/gentrome.fa.gz" >&2
    exit 1
  fi
  echo "[get_refs] decoys.txt: $(wc -l < "${decoys_txt}") sequences; gentrome.fa.gz ready" >&2
fi

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
  if [[ "${decoy}" == "true" ]]; then echo "- GENOME URL (decoy): ${genome_url}"; fi
  echo ""
  echo "## SHA256"
  echo "- $(basename "${gtf_gz}"): $(shasum -a 256 "${gtf_gz}" | awk '{print $1}')"
  echo "- $(basename "${fa_gz}"):  $(shasum -a 256 "${fa_gz}" | awk '{print $1}')"
  echo "- $(basename "${gtf_plain}"): $(shasum -a 256 "${gtf_plain}" | awk '{print $1}')"
  echo "- $(basename "${fa_plain}"):  $(shasum -a 256 "${fa_plain}" | awk '{print $1}')"
  if [[ "${decoy}" == "true" ]]; then
    echo "- $(basename "${genome_gz}"): $(shasum -a 256 "${genome_gz}" | awk '{print $1}')"
    echo "- decoys.txt: $(shasum -a 256 "${decoys_txt}" | awk '{print $1}')"
    echo "- gentrome.fa.gz: $(shasum -a 256 "${gentrome_gz}" | awk '{print $1}')"
  fi
} > "${readme}"

echo "[get_refs] Done. Files in: ${REFS_DIR}" >&2
echo "[get_refs] Contents:" >&2
ls -lh "${GTF_DIR}" "${FA_DIR}" >&2 || true
