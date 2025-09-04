#!/usr/bin/env bash
set -euo pipefail

# Usage: ./salmon_pipeline.sh [all|qc|trim|salmon]
#   all    = raw FastQC → trim → trimmed FastQC → MultiQC → Salmon index & quant
#   qc     = raw FastQC only
#   trim   = trim → trimmed FastQC → MultiQC → Salmon index & quant
#   salmon = Salmon index (if missing) → quant

START_STEP=${1:-all}
if [[ "$START_STEP" != "all" && "$START_STEP" != "qc" && "$START_STEP" != "trim" && "$START_STEP" != "salmon" ]]; then
  echo "Usage: $0 [all|qc|trim|salmon]"
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
    conda create -y -n "${ENV_NAME}" fastqc multiqc salmon trimmomatic -c bioconda -c conda-forge
  fi
  echo "Activating Conda env: ${ENV_NAME}"
  conda activate "${ENV_NAME}"
fi

# 3) Ensure output dirs exist
mkdir -p \
  "${FASTQC_RAW_DIR}" \
  "${TRIMMED_DIR}" \
  "${FASTQC_TRIM_DIR}" \
  "${MULTIQC_DIR}" \
  "${SALMON_OUT_DIR}" \
  "${LOGS_DIR}"

# 4) Checks
if [[ "$START_STEP" != "qc" && ! -d "$RAW_DIR" ]]; then
  echo "ERROR: Input FASTQ folder not found: $RAW_DIR" >&2
  exit 1
fi

echo "[05%] Setup complete"

# 5) FastQC on raw reads
if [[ "$START_STEP" == "all" || "$START_STEP" == "qc" ]]; then
  echo "[15%] FastQC (raw)"
  fastqc -t "${THREADS}" -o "${FASTQC_RAW_DIR}" "${RAW_DIR}"/*.fastq.gz
fi

# 6) Trim + FastQC (trimmed) + MultiQC
if [[ "$START_STEP" == "all" || "$START_STEP" == "trim" ]]; then
  echo "[35%] Trimmomatic"
  ADAPTER_FILE=( "${CONDA_PREFIX:-}"/share/trimmomatic-*/adapters/TruSeq3-PE.fa )
  if [[ ! -f "${ADAPTER_FILE[0]}" ]]; then
    echo "WARNING: Adapter file not found in conda env; using default name TruSeq3-PE.fa if available in CWD." >&2
  fi

  for R1 in "${RAW_DIR}"/*_R1_001.fastq.gz; do
    [[ -e "$R1" ]] || { echo "No FASTQs found in ${RAW_DIR}" >&2; break; }
    SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
    R2="${RAW_DIR}/${SAMPLE}_R2_001.fastq.gz"
    if [[ ! -f "$R2" ]]; then
      echo "WARNING: Mate not found for $SAMPLE; skipping." | tee -a "${LOGS_DIR}/pipeline.log"
      continue
    fi

    trimmomatic PE -threads "${THREADS}" \
      "$R1" "$R2" \
      "${TRIMMED_DIR}/${SAMPLE}_R1_trimmed.fastq.gz" \
      "${TRIMMED_DIR}/${SAMPLE}_R1_unpaired.fastq.gz" \
      "${TRIMMED_DIR}/${SAMPLE}_R2_trimmed.fastq.gz" \
      "${TRIMMED_DIR}/${SAMPLE}_R2_unpaired.fastq.gz" \
      ILLUMINACLIP:"${ADAPTER_FILE[0]:-TruSeq3-PE.fa}":2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
      2> "${LOGS_DIR}/trimmomatic_${SAMPLE}.log"
  done

  echo "[50%] FastQC (trimmed)"
  fastqc -t "${THREADS}" -o "${FASTQC_TRIM_DIR}" "${TRIMMED_DIR}"/*_trimmed.fastq.gz

  echo "[60%] MultiQC summary"
  multiqc "${OUT_DIR}" "${LOGS_DIR}" -o "${MULTIQC_DIR}"
fi

# 7) Salmon index & quantification
if [[ "$START_STEP" == "all" || "$START_STEP" == "salmon" ]]; then
  # Find FASTA (cdna) in data/references/fa
  FA_GZ=( "${FA_DIR}"/*.fa.gz )
  FA=( "${FA_DIR}"/*.fa )
  REF_FASTA=""
  if [[ -f "${FA_GZ[0]:-}" ]]; then REF_FASTA="${FA_GZ[0]}"; fi
  if [[ -z "$REF_FASTA" && -f "${FA[0]:-}" ]]; then REF_FASTA="${FA[0]}"; fi
  if [[ -z "$REF_FASTA" ]]; then
    echo "ERROR: No FASTA found in ${FA_DIR}. Use scripts/get_refs.sh first." >&2
    exit 1
  fi

  if [[ ! -d "${SALMON_INDEX}" ]]; then
    echo "[80%] Salmon index"
    salmon index -t "${REF_FASTA}" -i "${SALMON_INDEX}" -p "${THREADS}"
  fi

  echo "[90%] Salmon quant"
  for R1 in "${TRIMMED_DIR}"/*_R1_trimmed.fastq.gz; do
    [[ -e "$R1" ]] || { echo "No trimmed reads found in ${TRIMMED_DIR}" >&2; break; }
    SAMPLE=$(basename "$R1" _R1_trimmed.fastq.gz)
    R2="${TRIMMED_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
    salmon quant \
      -i "${SALMON_INDEX}" -l A \
      -1 "$R1" -2 "$R2" \
      -p "${THREADS}" \
      --gcBias --validateMappings \
      -o "${SALMON_OUT_DIR}/${SAMPLE}"
  done
fi

echo "[100%] Done"
echo "  • FastQC (raw):     ${FASTQC_RAW_DIR}"
echo "  • Trimmed reads:    ${TRIMMED_DIR}"
echo "  • FastQC (trimmed): ${FASTQC_TRIM_DIR}"
echo "  • MultiQC:          ${MULTIQC_DIR}"
echo "  • Salmon outputs:   ${SALMON_OUT_DIR}/<sample>/"
