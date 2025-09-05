#!/usr/bin/env bash
set -euo pipefail

# rename_fastqs.sh — Standardize FASTQ filenames for the pipeline
# Converts SRA-style pairs (e.g., SAMPLE_1.fastq.gz / SAMPLE_2.fastq.gz or .fq.gz)
# into the pipeline's expected pattern: SAMPLE_R1_001.fastq.gz / SAMPLE_R2_001.fastq.gz
#
# Usage examples:
#   bash scripts/rename_fastqs.sh                # rename files in data/fastq
#   bash scripts/rename_fastqs.sh --dry-run      # preview changes only
#   bash scripts/rename_fastqs.sh --dir path/to/fastqs
#
# Notes:
# - Only renames files that end with _1.fastq.gz, _2.fastq.gz, _1.fq.gz, or _2.fq.gz
# - Files already matching *_R1_001.fastq.gz or *_R2_001.fastq.gz are left unchanged
# - .fq.gz inputs are normalized to .fastq.gz outputs

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_DIR="${ROOT_DIR}/data/fastq"

target_dir="${DEFAULT_DIR}"
dry_run=0

usage() {
  cat >&2 <<EOF
Usage: $0 [--dir PATH] [--dry-run]

Renames SRA-style paired FASTQs to pipeline naming:
  SAMPLE_1.fastq.gz  -> SAMPLE_R1_001.fastq.gz
  SAMPLE_2.fastq.gz  -> SAMPLE_R2_001.fastq.gz
  SAMPLE_1.fq.gz     -> SAMPLE_R1_001.fastq.gz
  SAMPLE_2.fq.gz     -> SAMPLE_R2_001.fastq.gz

Options:
  --dir PATH   Directory containing FASTQs (default: ${DEFAULT_DIR})
  --dry-run    Show planned changes without renaming
  -h, --help   Show this message
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dir) target_dir="$2"; shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ ! -d "${target_dir}" ]]; then
  echo "[rename_fastqs][ERROR] Directory not found: ${target_dir}" >&2
  exit 1
fi

shopt -s nullglob

declare -a candidates
candidates=(
  "${target_dir}"/*_1.fastq.gz
  "${target_dir}"/*_2.fastq.gz
  "${target_dir}"/*_1.fq.gz
  "${target_dir}"/*_2.fq.gz
)

to_change=0
conflicts=0

rename_file() {
  local src="$1"
  local fname base read dest
  fname=$(basename "$src")

  # Skip if already in desired format
  if [[ "$fname" =~ _R[12]_001\.fastq\.gz$ ]]; then
    return 0
  fi

  if [[ "$fname" =~ ^(.+)_1\.(fastq|fq)\.gz$ ]]; then
    base="${BASH_REMATCH[1]}"; read=1
  elif [[ "$fname" =~ ^(.+)_2\.(fastq|fq)\.gz$ ]]; then
    base="${BASH_REMATCH[1]}"; read=2
  else
    echo "[rename_fastqs][SKIP] $fname (does not match *_1/*_2 patterns)" >&2
    return 0
  fi

  dest="${target_dir}/${base}_R${read}_001.fastq.gz"
  if [[ -e "$dest" && "$src" != "$dest" ]]; then
    echo "[rename_fastqs][ERROR] Destination exists: $(basename "$dest") (from $fname)" >&2
    conflicts=$((conflicts+1))
    return 1
  fi

  if (( dry_run )); then
    echo "DRY-RUN: $fname -> $(basename "$dest")"
  else
    mv -v "$src" "$dest"
  fi
  to_change=$((to_change+1))
}

echo "[rename_fastqs] Scanning: ${target_dir}"

found_any=0
for f in "${candidates[@]}"; do
  [[ -e "$f" ]] || continue
  found_any=1
  rename_file "$f" || true
done

if (( ! found_any )); then
  echo "[rename_fastqs] No SRA-style files found to rename in ${target_dir}." >&2
  exit 0
fi

if (( dry_run )); then
  echo "[rename_fastqs] Dry run complete. Planned renames: ${to_change}, conflicts: ${conflicts}" >&2
else
  echo "[rename_fastqs] Done. Renamed: ${to_change}, conflicts: ${conflicts}" >&2
fi

if (( conflicts > 0 )); then
  echo "[rename_fastqs][WARN] Some targets already existed. Resolve conflicts and rerun if needed." >&2
  exit 2
fi

exit 0

