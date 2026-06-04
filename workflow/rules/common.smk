"""Shared config, sample discovery, helpers, and component target rules."""

import os
import shutil


# ---------------------------------------------------------------------------
# Config-derived paths
# ---------------------------------------------------------------------------
FASTQ_DIR    = config["fastq_dir"]
OUT          = config["out_dir"]
THREADS      = int(config["threads"])
SALMON_INDEX = config["salmon_index"]

MULTIQC_REPORT = os.path.join(OUT, "multiqc", "multiqc_report.html")
QC_SUMMARY     = os.path.join(OUT, "salmon", "qc_summary.txt")


# ---------------------------------------------------------------------------
# Sample discovery
# Paired FASTQs named <sample>_R1_001.fastq.gz / <sample>_R2_001.fastq.gz in
# fastq_dir (same convention as salmon_pipeline.sh). Standardize SRA-style names
# first with scripts/rename_fastqs.sh (the fetch_test_data rule does this for you).
# ---------------------------------------------------------------------------
SAMPLES = sorted(set(glob_wildcards(os.path.join(FASTQ_DIR, "{sample}_R1_001.fastq.gz")).sample))


wildcard_constraints:
    sample=r"[^/]+",
    r=r"[12]",


# ---------------------------------------------------------------------------
# Target file lists
# ---------------------------------------------------------------------------
def trimmed_files():
    return expand(
        os.path.join(OUT, "trimmed", "{sample}_R{r}_trimmed.fastq.gz"),
        sample=SAMPLES, r=["1", "2"],
    )


def fastqc_files():
    return expand(
        os.path.join(OUT, "fastqc_trimmed", "{sample}_R{r}_trimmed_fastqc.html"),
        sample=SAMPLES, r=["1", "2"],
    )


def quant_files():
    return expand(os.path.join(OUT, "salmon", "{sample}", "quant.sf"), sample=SAMPLES)


def deseq2_target():
    pn = config["deseq2"]["project_name"]
    return os.path.join(OUT, "deseq2", f"{pn}_DESeq2_full_results.csv")


def all_targets():
    """Default `rule all` targets: MultiQC + Salmon QC gate (+ DESeq2 unless disabled)."""
    targets = [MULTIQC_REPORT, QC_SUMMARY]
    if config["deseq2"].get("run", True):
        targets.append(deseq2_target())
    return targets


# ---------------------------------------------------------------------------
# Rule helpers
# ---------------------------------------------------------------------------
def salmon_index_inputs(wildcards):
    """Decoy-aware -> gentrome + decoys; otherwise transcriptome only. All ancient()
    so a refreshed reference timestamp can never trigger an index rebuild."""
    if config.get("decoy_aware", True):
        return {
            "ref": ancient(config["gentrome"]),
            "decoys": ancient(config["decoys"]),
        }
    return {"ref": ancient(config["transcriptome"])}


def deseq2_optional_args(wildcards):
    """Build the optional run_deseq2.R flags from config, omitting null entries."""
    d = config["deseq2"]
    parts = []
    if d.get("design"):
        parts.append("--design '{}'".format(d["design"]))
    if d.get("test"):
        parts.append("--test {}".format(d["test"]))
    if d.get("reduced"):
        parts.append("--reduced '{}'".format(d["reduced"]))
    if d.get("ref_level"):
        parts.append("--ref_level {}".format(d["ref_level"]))
    if d.get("contrast"):
        parts.append("--contrast {}".format(d["contrast"]))
    return " ".join(parts)


# ---------------------------------------------------------------------------
# Component target rules — run any stage on its own
#   snakemake trim_all | fastqc_all | quant_all | qc_check | multiqc | deseq2 | index | clean
# ---------------------------------------------------------------------------
localrules: trim_all, fastqc_all, quant_all, qc_check, index, clean


rule trim_all:
    input:
        trimmed_files(),


rule fastqc_all:
    input:
        fastqc_files(),


rule quant_all:
    input:
        quant_files(),


rule qc_check:
    input:
        QC_SUMMARY,


rule index:
    input:
        SALMON_INDEX,


rule clean:
    # Removes ONLY generated pipeline outputs (out/ and logs/). Deliberately NEVER
    # touches: data/references/, data/fastq/ (raw data), salmon_index/ or
    # resources/salmon_index_* (the indices), or .snakemake/ (metadata + conda
    # cache + the prebuilt-index registration). Fixed allow-list, no globbing.
    run:
        for path in [OUT, "logs"]:
            if os.path.isdir(path):
                shutil.rmtree(path)
                print(f"[clean] removed {path}/")
            else:
                print(f"[clean] {path}/ absent, nothing to do")
