#!/usr/bin/env bash
set -euo pipefail

# fetch_test_fastqs.sh — Fetch FASTQs for a GEO series into data/fastq/
# Default dataset: GSE52778

# Cleanup temporary files on exit
cleanup() {
    find /tmp -maxdepth 1 -name "fetch_test_fastqs.*" -type f -exec rm -f {} +
}
trap cleanup EXIT INT TERM
#
# Usage examples:
#   bash scripts/fetch_test_fastqs.sh
#   bash scripts/fetch_test_fastqs.sh --geo GSE52778 --method ena
#   bash scripts/fetch_test_fastqs.sh --runs SRR123,SRR456 --out data/fastq
#   THREADS=8 bash scripts/fetch_test_fastqs.sh
#
# Notes:
# - Default method is 'sra-tools' (uses SRA Toolkit's fasterq-dump).
#   You can force --method ena to fetch from ENA instead.
# - Requires network access. ENA downloads resume with curl -C -.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
FASTQ_DIR="${ROOT_DIR}/data/fastq"  # Base fastq directory

geo="GSE52778"
runs_csv=""
out_dir=""  # Will be set to FASTQ_DIR/subdir if --out specifies a subdir
method="sra-tools"   # ena|sra-tools|auto (default: sra-tools)
threads="${THREADS:-4}"
parallel_downloads="${PARALLEL:-4}"  # Number of parallel downloads, defaults to 4

usage() {
  echo "Usage: $0 [--geo GSE52778] [--runs SRR1,SRR2,...] [--out SUBDIR] [--method ena|sra-tools|auto] [--threads N] [--parallel N]" >&2
  echo "Note: Output files will always be placed under ${FASTQ_DIR}. The --out parameter specifies an optional subdirectory." >&2
  echo "      Use --parallel N to set number of concurrent downloads (default: 4, or set PARALLEL env var)" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --geo) geo="$2"; shift 2 ;;
    --runs) runs_csv="$2"; shift 2 ;;
    --out) 
      # If absolute path provided, ensure it's under data/fastq
      if [[ "$2" = /* ]]; then
        if ! [[ "$2" =~ ^"${FASTQ_DIR}"(/|$) ]]; then
          err "Output directory must be under ${FASTQ_DIR}"
          exit 1
        fi
        out_dir="$2"
      else
        # Relative path - append to FASTQ_DIR
        out_dir="${FASTQ_DIR}/${2#/}"
      fi
      shift 2 ;;
    --method) method="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --parallel) parallel_downloads="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

# If no output directory specified, use base FASTQ_DIR
if [[ -z "${out_dir}" ]]; then
  out_dir="${FASTQ_DIR}"
fi

# Validate output directory
if ! mkdir -p "${out_dir}" 2>/dev/null; then
    err "Cannot create output directory: ${out_dir}"
    exit 1
fi
if ! [[ -w "${out_dir}" ]]; then
    err "Output directory is not writable: ${out_dir}"
    exit 1
fi

# Check available disk space (need at least 10GB free)
if command -v df >/dev/null 2>&1; then
    free_space=$(df -k "${out_dir}" | awk 'NR==2 {print $4}')
    if [[ ${free_space} -lt 10485760 ]]; then  # 10GB in KB
        err "Insufficient disk space. Need at least 10GB free in ${out_dir}"
        exit 1
    fi
fi

log() { echo "[fetch_test] $*" >&2; }
err() { echo "[fetch_test][ERROR] $*" >&2; }

resolve_runs_from_geo() {
  local acc="$1"
  # Validate GEO accession format
  if ! [[ $acc =~ ^GSE[0-9]+$ ]]; then
    err "Invalid GEO accession format: ${acc} (should be GSExxxxx)"
    return 1
  fi

  log "Resolving SRR runs from GEO: ${acc}"

  # Step 1: Get SRP accession from GEO
  local geo_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${acc}&targ=self&form=text&view=quick"
  local tmp_geo tmp_ena
  tmp_geo="$(mktemp)"

  if ! curl -fsSL -o "$tmp_geo" "$geo_url"; then
    rm -f "$tmp_geo"
    err "Failed to retrieve GEO record for ${acc}."
    return 1
  fi

  # Extract SRP accession from GEO record
  local srp
  srp=$(grep -i "Series_relation.*SRA.*SRP" "$tmp_geo" | grep -oE 'SRP[0-9]+' | head -n1)
  rm -f "$tmp_geo"

  if [[ -z "$srp" ]]; then
    err "No SRA project (SRP) found for ${acc}. This GEO series may not have sequencing data."
    return 1
  fi

  log "Found SRA project: ${srp}"

  # Step 2: Get run accessions from ENA using the SRP
  local ena_url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srp}&result=read_run&fields=run_accession&format=tsv"
  tmp_ena="$(mktemp)"

  if ! curl -fsSL -o "$tmp_ena" "$ena_url"; then
    rm -f "$tmp_ena"
    err "Failed to retrieve run list from ENA for ${srp}."
    return 1
  fi

  # Extract run accessions (skip header)
  local runs
  runs=$(tail -n +2 "$tmp_ena" | tr -d '\r' | sort -u)
  rm -f "$tmp_ena"

  if [[ -z "$runs" ]]; then
    err "No runs found for ${srp}."
    return 1
  fi

  echo "$runs"
}

download_run_via_ena() {
  local run="$1"; local out="$2"
  local api="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv&download=false"
  local tmp
  tmp="$(mktemp)"
  if ! curl -fsSL -o "$tmp" "$api"; then
    rm -f "$tmp"
    return 1
  fi
  # Parse TSV: header then one row
  local line
  line="$(tail -n +2 "$tmp" | head -n1)"
  rm -f "$tmp"
  if [[ -z "$line" ]]; then
    return 1
  fi
  IFS=$'\t' read -r run_acc fastq_ftp fastq_md5 <<< "$line"
  if [[ -z "$fastq_ftp" || "$fastq_ftp" == "null" ]]; then
    return 1
  fi

  IFS=';' read -r -a files <<< "$fastq_ftp"
  IFS=';' read -r -a md5s <<< "$fastq_md5"

  for i in "${!files[@]}"; do
    local f="${files[$i]}"
    local url="https://${f}"
    local fname
    fname=$(basename "$f")
    local dest="${out}/${fname}"
    log "Downloading ${run}: ${url} -> ${dest}"
    curl -fL -C - -o "$dest" "$url"
    if [[ ! -s "$dest" ]]; then
      err "Downloaded file is empty: $dest"
      return 1
    fi
    # Optional MD5 verification where tool exists and MD5 provided
    local expected="${md5s[$i]:-}"
    if [[ -n "$expected" && "$expected" != "null" ]]; then
      local got=""
      if command -v md5sum >/dev/null 2>&1; then
        got=$(md5sum "$dest" | awk '{print $1}')
      elif command -v md5 >/dev/null 2>&1; then
        got=$(md5 -q "$dest")
      fi
      if [[ -n "$got" && "$got" != "$expected" ]]; then
        err "MD5 mismatch for $fname (got $got expected $expected)"
        return 1
      fi
    fi
  done
}

download_run_via_sra_tools() {
  local run="$1"; local out="$2"; local thr="$3"
  if ! command -v fasterq-dump >/dev/null 2>&1; then
    err "SRA Toolkit not found. Install with: mamba install -c bioconda sra-tools"
    return 1
  fi
  log "Using SRA Toolkit fasterq-dump for ${run} (threads=${thr})"
  fasterq-dump "${run}" -e "${thr}" -O "$out"
  # Compress outputs for consistency (use pigz if available)
  local compressor="gzip -f"
  if command -v pigz >/dev/null 2>&1; then
    compressor="pigz -f -p ${thr}"
  fi
  find "$out" -maxdepth 1 -type f -name "${run}_*.fastq" -print0 | xargs -0 -I{} sh -c "${compressor} \"{}\""
}

# 1) Resolve SRR list
runs_list=()
if [[ -n "$runs_csv" ]]; then
  IFS=',' read -r -a runs_list <<< "$runs_csv"
else
  while IFS= read -r run; do
    runs_list+=("$run")
  done < <(resolve_runs_from_geo "$geo" || true)
fi

if [[ ${#runs_list[@]} -eq 0 ]]; then
  err "No runs resolved. Pass --runs SRR1,SRR2 or check GEO accession."
  exit 1
fi

log "Runs to fetch (${#runs_list[@]}): ${runs_list[*]}"

# Function to process a single run (for parallel execution)
process_run() {
    local r="$1"
    local out="$2"
    local meth="$3"
    local thr="$4"
    
    if [[ "$meth" == "ena" || "$meth" == "auto" ]]; then
        if download_run_via_ena "$r" "$out"; then
            echo "[fetch_test] Successfully downloaded ${r} via ENA"
            return 0
        elif [[ "$meth" == "ena" ]]; then
            echo "[fetch_test][ERROR] ENA download failed for ${r}"
            return 1
        fi
    fi
    if [[ "$meth" == "sra-tools" || "$meth" == "auto" ]]; then
        if download_run_via_sra_tools "$r" "$out" "$thr"; then
            echo "[fetch_test] Successfully downloaded ${r} via SRA Tools"
            return 0
        else
            echo "[fetch_test][ERROR] SRA Toolkit download failed for ${r}"
            return 1
        fi
    fi
    return 1
}

# 2) Download runs in parallel
if command -v parallel >/dev/null 2>&1; then
    log "Using GNU parallel for downloads (jobs: ${parallel_downloads})"
    
    # Export functions and variables needed by parallel
    export -f process_run download_run_via_ena download_run_via_sra_tools log err
    export method threads out_dir
    
    # Create a temporary file to store results
    results_file="$(mktemp /tmp/fetch_test_fastqs.results.XXXXXX)"
    
    # Run downloads in parallel
    printf "%s\n" "${runs_list[@]}" | \
        parallel --will-cite -j "${parallel_downloads}" \
        "process_run {} '${out_dir}' '${method}' '${threads}' || echo FAIL_{}" > "${results_file}"
    
    # Count failures
    fail_count=$(grep -c '^FAIL_' "${results_file}" || true)
    rm -f "${results_file}"
else
    log "GNU parallel not found, running downloads sequentially"
    fail_count=0
    for r in "${runs_list[@]}"; do
        if ! process_run "$r" "$out_dir" "$method" "$threads"; then
            fail_count=$((fail_count+1))
        fi
    done
fi

log "Done. FASTQs in: ${out_dir}"
if [[ $fail_count -gt 0 ]]; then
  err "There were $fail_count failures. See logs above."
  exit 2
fi
