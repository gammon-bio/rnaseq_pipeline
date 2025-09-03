# Salmon + DESeq2 Pipeline (with your Rmd logic)

This repository packages your working Salmon quantification and DESeq2 downstream analysis into a reproducible, GitHub-ready layout with minimal, surgical changes. Your original analysis logic (volcano, enrichment, etc.) remains intact inside `tximport_deseq2.rmd`; the only additions are parameters for input/output paths and thresholds so it can be run headlessly.

## Quickstart

- Create environments (CLI tools and R):
  - mamba env create -f environment.yml
  - mamba env create -f environment-r.yml

- Activate CLI env and fetch references:
  - conda activate rnaseq
  - bash scripts/get_refs.sh --species human --build GRCh38
    - For mouse: `--species mouse --build GRCm39`
    - To pin Ensembl release: `--release 110`
    - Or pass explicit URLs: `--gtf_url ... --fasta_url ...`

- Run the Salmon pipeline (raw QC → trim → trimmed QC → MultiQC → quant):
  - bash salmon_pipeline.sh all
  - Inputs: `data/fastq/` with `*_R1_001.fastq.gz` and `*_R2_001.fastq.gz`
  - Outputs:
    - out/fastqc_raw/
    - out/trimmed/
    - out/fastqc_trimmed/
    - out/multiqc/
    - out/salmon/<sample>/ (quant.sf)
    - logs/

- Activate R env and run DESeq2 wrapper (renders your Rmd headlessly):
  - conda activate rnaseq-r
  - Rscript scripts/run_deseq2.R \
      --quant_dir out/salmon \
      --gtf data/references/gtf/<your>.gtf \
      --sample_table examples/sample_table.csv \
      --group_col condition \
      --project_name CU25 \
      --padj_thresh 0.05 --lfc_thresh 0.5
  - Outputs (written by your Rmd under `out/deseq2/`):
    - DE results CSV (e.g., `DESeq2_full_results_*.csv`)
    - Volcano PDF (e.g., `volcano_plot_*.pdf`)
    - PCA PDF (e.g., `PCA_plot_*.pdf`)
    - VST normalized counts CSV (e.g., `vst_norm_counts_*.csv`)
    - Up/Down enrichment Excel files if your Rmd writes them

## References retrieval (scripts/get_refs.sh)

- Flags:
  - `--species` human|mouse (default: human)
  - `--build` GRCh38|GRCm39 (default depends on species)
  - `--release` <n>|current (default: current)
  - `--gtf_url` and `--fasta_url` to override URLs directly
- Behavior:
  - Creates `refs/{gtf,fa}/`
  - Downloads via `curl -L -C -` (resume)
  - Decompresses `.gz` to `.gtf`/`.fa`
  - Writes `refs/README.md` with exact URLs and SHA256 checksums
- Salmon references use Ensembl cDNA FASTA (best practice for transcript-level quantification).

## DESeq2 runner (scripts/run_deseq2.R)

- Modes:
  - `--quant_dir + --gtf` → Rmd builds `tx2gene` from GTF and imports Salmon (default)
  - `--quant_dir + --tx2gene` → wrapper computes `tximport` and passes `--tximport_rds` to the Rmd
  - `--tximport_rds` → use a precomputed `tximport` object directly
- Common flags:
  - `--sample_table` CSV with sample metadata (row names = sample IDs)
  - `--group_col` design column (default: `condition`)
  - `--padj_thresh`, `--lfc_thresh` forwarded to your volcano/summary logic
  - `--out_dir` output directory (default: `out/deseq2`)

## Notes on practices and reproducibility

- This pipeline aligns with Bioconductor’s RNA-seq gene-level workflow guidance for design, dispersion, and multiple testing.
  - Bioconductor workflow: https://bioconductor.org/packages/release/workflows/html/rnaseqGene.html
- Please cite the tools used:
  - Salmon: Patro et al., Nature Methods 2017
  - DESeq2: Love et al., Genome Biology 2014
  - tximport: Soneson et al., F1000Research 2015
  - MultiQC: Ewels et al., Bioinformatics 2016
  - Ensembl/biomaRt for annotation

## Repository layout

- salmon_pipeline/
  - README.md, LICENSE, .gitignore
  - environment.yml (CLI), environment-r.yml (R)
  - salmon_pipeline.sh
  - tximport_deseq2.rmd
  - scripts/
    - get_refs.sh
    - run_deseq2.R
  - data/
    - fastq/ (place raw FASTQs here)
    - references/ (created by get_refs.sh)
  - out/ (pipeline outputs)
  - examples/
    - sample_table.csv

## Backward compatibility

- `salmon_pipeline.sh` keeps your original entrypoint and step flag (`all|qc|trim|salmon`). It now logs progress with percentages and writes outputs to `out/` subfolders. MultiQC is pointed at FastQC outputs and Trimmomatic logs under `logs/` so it picks them up.
- `tximport_deseq2.rmd` is unchanged in analysis logic; only parameters were added for file paths and thresholds so your exact volcano and enrichment code is preserved.
