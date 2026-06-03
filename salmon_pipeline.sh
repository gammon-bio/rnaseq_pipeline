#!/usr/bin/env bash
set -euo pipefail

# Usage: ./salmon_pipeline.sh [all|trim|salmon]
#   all    = fastp trim → trimmed FastQC → MultiQC → Salmon index & quant
#   trim   = fastp trim → trimmed FastQC → MultiQC → Salmon index & quant
#   salmon = Salmon index (if missing) → quant

START_STEP=${1:-all}
if [[ "$START_STEP" != "all" && "$START_STEP" != "trim" && "$START_STEP" != "salmon" ]]; then
  echo "Usage: $0 [all|trim|salmon]"
  exit 1
fi

# 1) Project directories (script‑relative)
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="${PROJECT_DIR}/data/fastq"                  # input FASTQs
OUT_DIR="${PROJECT_DIR}/out"
FASTQC_RAW_DIR="${OUT_DIR}/fastqc_raw"
TRIMMED_DIR="${OUT_DIR}/trimmed"
FASTQC_TRIM_DIR="${OUT_DIR}/fastqc_trimmed"
MULTIQC_DIR="${OUT_DIR}/multiqc"
SALMON_OUT_DIR="${OUT_DIR}/salmon"
LOGS_DIR="${PROJECT_DIR}/logs"

# References
REFS_DIR="${PROJECT_DIR}/data/references"
FA_DIR="${REFS_DIR}/fa"
GTF_DIR="${REFS_DIR}/gtf"
SALMON_INDEX="${PROJECT_DIR}/salmon_index"

THREADS=${THREADS:-8}
ENV_NAME=rnaseq

# 2) Conda activate if available (backward‑compatible)
if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
  # shellcheck source=/dev/null
  source "${CONDA_BASE}/etc/profile.d/conda.sh"
  if ! conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    echo "Creating Conda env: ${ENV_NAME}"
    conda create -y -n "${ENV_NAME}" fastqc multiqc salmon fastp -c conda-forge -c bioconda
  fi
  echo "Activating Conda env: ${ENV_NAME}"
  conda activate "${ENV_NAME}"
fi

# 3) Ensure output dirs exist
mkdir -p \
  "${TRIMMED_DIR}" \
  "${FASTQC_TRIM_DIR}" \
  "${MULTIQC_DIR}" \
  "${SALMON_OUT_DIR}" \
  "${LOGS_DIR}"

# 4) Validate input directories based on mode
if [[ "$START_STEP" == "all" || "$START_STEP" == "trim" ]]; then
  if [[ ! -d "$RAW_DIR" ]]; then
    echo "ERROR: Input FASTQ folder not found: $RAW_DIR" >&2
    exit 1
  fi
elif [[ "$START_STEP" == "salmon" ]]; then
  if [[ ! -d "$TRIMMED_DIR" ]]; then
    echo "ERROR: Trimmed reads folder not found: $TRIMMED_DIR" >&2
    echo "Run 'bash salmon_pipeline.sh trim' first, or place trimmed reads in ${TRIMMED_DIR}" >&2
    exit 1
  fi
fi

echo "[10%] Setup complete"

# 5) Trim with fastp + FastQC (trimmed)
if [[ "$START_STEP" == "all" || "$START_STEP" == "trim" ]]; then
  echo "[25%] fastp trimming and filtering"

  for R1 in "${RAW_DIR}"/*_R1_001.fastq.gz; do
    [[ -e "$R1" ]] || { echo "No FASTQs found in ${RAW_DIR}" >&2; break; }
    SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
    R2="${RAW_DIR}/${SAMPLE}_R2_001.fastq.gz"
    if [[ ! -f "$R2" ]]; then
      echo "WARNING: Mate not found for $SAMPLE; skipping." | tee -a "${LOGS_DIR}/pipeline.log"
      continue
    fi

    fastp \
      -i "$R1" -I "$R2" \
      -o "${TRIMMED_DIR}/${SAMPLE}_R1_trimmed.fastq.gz" \
      -O "${TRIMMED_DIR}/${SAMPLE}_R2_trimmed.fastq.gz" \
      --trim_poly_g \
      --qualified_quality_phred 20 \
      --length_required 36 \
      --thread "${THREADS}" \
      --json "${LOGS_DIR}/${SAMPLE}.fastp.json" \
      --html "${LOGS_DIR}/${SAMPLE}.fastp.html" \
      2>&1 | tee "${LOGS_DIR}/fastp_${SAMPLE}.log"
  done

  echo "[45%] FastQC (trimmed)"
  TRIMMED_FASTQS=( "${TRIMMED_DIR}"/*_trimmed.fastq.gz )
  if [[ ! -f "${TRIMMED_FASTQS[0]:-}" ]]; then
    echo "ERROR: No trimmed .fastq.gz files found in ${TRIMMED_DIR}" >&2
    exit 1
  fi
  fastqc -t "${THREADS}" -o "${FASTQC_TRIM_DIR}" "${TRIMMED_FASTQS[@]}"
fi

# 6) Salmon index & quantification
if [[ "$START_STEP" == "all" || "$START_STEP" == "salmon" ]]; then
  GENTROME="${FA_DIR}/gentrome.fa.gz"
  DECOYS="${FA_DIR}/decoys.txt"

  # Decide whether the Salmon index needs to be (re)built
  NEED_INDEX=false
  if [[ ! -d "${SALMON_INDEX}" ]]; then
    NEED_INDEX=true
  elif [[ -f "${GENTROME}" && -f "${DECOYS}" ]] && grep -q '"num_decoys": 0' "${SALMON_INDEX}/info.json" 2>/dev/null; then
    echo "Existing Salmon index has no decoys but decoy references are present; rebuilding."
    rm -rf "${SALMON_INDEX}"
    NEED_INDEX=true
  fi

  if [[ "${NEED_INDEX}" == "true" ]]; then
    if [[ -f "${GENTROME}" && -f "${DECOYS}" ]]; then
      echo "[65%] Salmon index (decoy-aware: gentrome + decoys)"
      salmon index -t "${GENTROME}" -d "${DECOYS}" -i "${SALMON_INDEX}" -k 31 -p "${THREADS}"
    else
      # Transcriptome-only fallback: prefer a cDNA FASTA, excluding genome/gentrome files
      REF_FASTA=""
      for f in "${FA_DIR}"/*cdna*.fa.gz "${FA_DIR}"/*.fa.gz "${FA_DIR}"/*.fa; do
        [[ -e "$f" ]] || continue
        case "$(basename "$f")" in
          gentrome.fa.gz|*dna.primary_assembly*|*dna.toplevel*|*dna_sm*|*dna_rm*) continue ;;
        esac
        REF_FASTA="$f"; break
      done
      if [[ -z "${REF_FASTA}" ]]; then
        echo "ERROR: No transcriptome FASTA found in ${FA_DIR}. Use scripts/get_refs.sh first." >&2
        exit 1
      fi
      echo "[65%] Salmon index (transcriptome-only: $(basename "${REF_FASTA}"))"
      salmon index -t "${REF_FASTA}" -i "${SALMON_INDEX}" -k 31 -p "${THREADS}"
    fi
  fi

  echo "[75%] Salmon quant"
  TRIMMED_R1=( "${TRIMMED_DIR}"/*_R1_trimmed.fastq.gz )
  if [[ ! -f "${TRIMMED_R1[0]:-}" ]]; then
    echo "ERROR: No trimmed R1 reads found in ${TRIMMED_DIR}" >&2
    exit 1
  fi

  for R1 in "${TRIMMED_R1[@]}"; do
    SAMPLE=$(basename "$R1" _R1_trimmed.fastq.gz)
    R2="${TRIMMED_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"

    if [[ ! -f "$R2" ]]; then
      echo "ERROR: Missing R2 mate for sample ${SAMPLE}: $R2" >&2
      exit 1
    fi

    salmon quant \
      -i "${SALMON_INDEX}" -l A \
      -1 "$R1" -2 "$R2" \
      -p "${THREADS}" \
      --gcBias --validateMappings \
      -o "${SALMON_OUT_DIR}/${SAMPLE}"
  done
fi

# 7) Final MultiQC (captures fastp, FastQC trimmed, and Salmon)
# Build MultiQC input list based on what was actually run
MULTIQC_INPUTS=()

if [[ -d "${FASTQC_TRIM_DIR}" && -n "$(ls -A "${FASTQC_TRIM_DIR}" 2>/dev/null)" ]]; then
  MULTIQC_INPUTS+=("${FASTQC_TRIM_DIR}")
fi

if [[ -d "${SALMON_OUT_DIR}" && -n "$(ls -A "${SALMON_OUT_DIR}" 2>/dev/null)" ]]; then
  MULTIQC_INPUTS+=("${SALMON_OUT_DIR}")
fi

if [[ -d "${LOGS_DIR}" && -n "$(ls -A "${LOGS_DIR}" 2>/dev/null)" ]]; then
  MULTIQC_INPUTS+=("${LOGS_DIR}")
fi

if [[ ${#MULTIQC_INPUTS[@]} -gt 0 ]]; then
  echo "[95%] MultiQC summary"
  multiqc "${MULTIQC_INPUTS[@]}" -o "${MULTIQC_DIR}" --force
else
  echo "[95%] Skipping MultiQC - no outputs to summarize"
fi

echo "[100%] Done"
echo "  • Trimmed reads:    ${TRIMMED_DIR}"
echo "  • FastQC (trimmed): ${FASTQC_TRIM_DIR}"
echo "  • MultiQC:          ${MULTIQC_DIR}"
echo "  • Salmon outputs:   ${SALMON_OUT_DIR}/<sample>/"
