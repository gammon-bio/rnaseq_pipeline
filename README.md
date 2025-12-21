# Salmon + DESeq2 Pipeline 

This repository combines RNAseq fastq file processing, Salmon and DESeq2 downstream analysis into a reproducible pipeline. Fastq files are processed with fastp, aligned/quantified with Salmon then analyzed with DESeq2. Final output includes annotated raw count files, DESeq2 statistics output table, differential experession visuals (PCA and volcano plot) and functional enrichment results. Contains multiple helper scripts for fetching fastqs from public datasets, fetching fa and gtf reference files and python based QC of salmon alignment.

## Quickstart

**1. Create environments (CLI tools and R):**

For latest versions:
```bash
mamba env create -f environment.yml
mamba env create -f environment-r.yml
```

For exact reproducibility:
```bash
mamba env create -f environment.lock.yml
mamba env create -f environment-r.lock.yml
```

**2. Fetch reference genome and transcriptome:**

Activate the CLI environment and download references:
```bash
conda activate rnaseq
bash scripts/get_refs.sh --species human --build GRCh38
```

For mouse data:
```bash
bash scripts/get_refs.sh --species mouse --build GRCm39
```

To pin a specific Ensembl release (e.g., release 110):
```bash
bash scripts/get_refs.sh --species human --build GRCh38 --release 110
```

Or pass explicit GTF/FASTA URLs:
```bash
bash scripts/get_refs.sh --gtf_url https://... --fasta_url https://...
```

**3. Prepare FASTQs and run Salmon pipeline:**

Place your paired-end reads in `data/fastq/` with naming: `<SAMPLE>_R1_001.fastq.gz` and `<SAMPLE>_R2_001.fastq.gz`

If your files follow SRA naming convention (`SAMPLE_1.fastq.gz` or `SAMPLE_1.fq.gz`), standardize them first:
```bash
bash scripts/rename_fastqs.sh
```

Run the full Salmon pipeline (fastp trim → trimmed QC → MultiQC → quant):
```bash
bash salmon_pipeline.sh all
```

To use multiple CPU cores (example: 8 cores):
```bash
THREADS=8 bash salmon_pipeline.sh all
```

**Performance note:** Pipeline uses fastp for trimming (~5x faster than Trimmomatic) with automatic adapter detection, poly-G tail trimming, and comprehensive HTML/JSON QC reports. Raw FastQC is skipped since fastp handles all QC issues (adapters, quality, poly-G); post-trim FastQC + Salmon metrics provide sufficient validation.

Pipeline outputs:
- `out/trimmed/` — Trimmed read pairs (fastp output)
- `out/fastqc_trimmed/` — FastQC reports on trimmed reads
- `out/multiqc/` — MultiQC summary report (includes fastp stats)
- `out/salmon/<sample>/` — Salmon quantification (quant.sf)
- `logs/` — Pipeline logs (includes fastp JSON/HTML reports)

### FASTQ Naming and Renaming

- Expected naming for paired-end reads: `<SAMPLE>_R1_001.fastq.gz` and `<SAMPLE>_R2_001.fastq.gz` in `data/fastq/`.
  - Note: fastp trimming only processes files matching `*_R1_001.fastq.gz` pattern, so failure to name properly will cause pipeline to skip samples
- If your files are SRA-style (e.g., `SAMPLE_1.fastq.gz` / `SAMPLE_2.fastq.gz` or `.fq.gz`), use the helper script to standardize names:
  - Preview changes: `bash scripts/rename_fastqs.sh --dry-run`
  - Apply changes: `bash scripts/rename_fastqs.sh`
  - Custom folder: `bash scripts/rename_fastqs.sh --dir path/to/fastqs`
  - The script converts `_1/_2` to `_R1_001/_R2_001` and normalizes `.fq.gz` to `.fastq.gz`.

### QC Check (Recommended)

After Salmon quantification, verify mapping rates before proceeding to differential expression:

```bash
python scripts/check_salmon_qc.py --salmon_dir out/salmon --min_mapping 0.6
```

Options:
- `--salmon_dir`: Path to Salmon output directory (default: `out/salmon`)
- `--min_mapping`: Minimum mapping rate threshold as decimal (default: `0.6` = 60%)
- `--out`: Optional path to save summary file (e.g., `--out out/salmon_qc_summary.txt`)

Samples below the mapping rate threshold will be flagged. **Remove failed samples from your `sample_table.csv` before running DESeq2.** The script exits with code 1 if any samples fail, allowing use in automated pipelines.

> **Note:** Low mapping rates (<50%) often indicate reference mismatch, adapter contamination, or sample quality issues. Check the MultiQC report for diagnostic details before excluding samples.

**4. Run differential expression analysis with DESeq2:**

Activate the R environment and run the DESeq2 wrapper:
```bash
conda activate rnaseq-r
Rscript scripts/run_deseq2.R \
  --quant_dir out/salmon \
  --gtf data/references/gtf/Homo_sapiens.GRCh38.115.gtf \
  --sample_table examples/sample_table.csv \
  --group_col condition \
  --project_name MyProject \
  --padj_thresh 0.05 \
  --lfc_thresh 0.5
```

For mouse data:
```bash
Rscript scripts/run_deseq2.R \
  --quant_dir out/salmon \
  --gtf data/references/gtf/Mus_musculus.GRCm39.115.gtf \
  --sample_table examples/sample_table.csv \
  --group_col condition \
  --project_name MyProject \
  --padj_thresh 0.05 \
  --lfc_thresh 0.5
```

**Important:** The `--gtf` flag requires the complete filename including extension (e.g., `Homo_sapiens.GRCh38.115.gtf`), not just the directory path.

### DESeq2 outputs

All outputs are written to `out/deseq2/` with your `--project_name` as prefix:

- `<PROJECT>_DESeq2_full_results.csv` — Full DE results with gene symbols (baseMean, log2FC, p-value, padj)
- `<PROJECT>_volcano_plot.pdf` — Volcano plot with significant genes labeled
- `<PROJECT>_PCA_plot.pdf` — PCA plot on VST-transformed counts
- `<PROJECT>_vst_norm_counts.csv` — VST-normalized counts with gene symbols
- `<PROJECT>_gene_counts_with_symbols.csv` — Raw counts with gene symbols
- `<PROJECT>_raw_gene_counts.csv` — Raw counts (Ensembl IDs only)
- `<PROJECT>_Upregulated_cachexia_Enrichment_results.xlsx` — Gene enrichment for upregulated genes
- `<PROJECT>_Downregulated_cachexia_Enrichment_results.xlsx` — Gene enrichment for downregulated genes

**Note:** Gene symbols are automatically detected based on your GTF file (human vs mouse) and merged into all count matrices and results tables.

## References retrieval (scripts/get_refs.sh)

- Flags:
  - `--species` human|mouse (default: human)
  - `--build` GRCh38|GRCm39 (default depends on species)
  - `--release` <n>|current (default: current)
  - `--gtf_flavor` plain|chr|chr_patch_hapl_scaff|abinitio|auto (default: auto; prefers plain)
  - `--gtf_url` and `--fasta_url` to override URLs directly
- Behavior:
  - Creates `data/references/{gtf,fa}/`
  - Downloads via `curl -L -C -` (resume)
  - Decompresses `.gz` to `.gtf`/`.fa`
  - Writes `data/references/README.md` with exact URLs and SHA256 checksums
- Salmon references use Ensembl cDNA FASTA (best practice for transcript-level quantification).

## Optional: Test Dataset (GSE52778)

GSE52778 is a human airway smooth muscle cell dataset with 4 samples (2 controls + 2 treatments). Perfect for testing the full pipeline end-to-end. Pipeline correctly identifies top DEGs (FKBP5, MAOA, KLF15) from PMID: 24926665

### Fetch test FASTQs

**Option 1: Use ENA instead of SRA Toolkit (faster, curl-based):**
```bash
bash scripts/fetch_test_fastqs.sh --runs SRR1039508,SRR1039509,SRR1039512,SRR1039513 --geo GSE52778 --method ena
```

**Option 2: Fetch by GEO accession (automatic run resolution, must install SRA-tools into conda env first):**
```bash
conda activate rnaseq
bash scripts/fetch_test_fastqs.sh --geo GSE52778
```

**Option 3: Fetch specific SRR runs directly:**
```bash
bash scripts/fetch_test_fastqs.sh --runs SRR1039508,SRR1039509,SRR1039512,SRR1039513
```

**Option 4: Use multiple threads with SRA Toolkit:**
```bash
bash scripts/fetch_test_fastqs.sh --geo GSE52778 --method sra-tools --threads 8
```

### Download methods explained

- `--method sra-tools` (default): Uses SRA Toolkit's `fasterq-dump`. Requires: `mamba install -c bioconda sra-tools`
- `--method ena`: Fetches from ENA using curl (faster, resumes with `-C -`). No additional tools needed.
- `--method auto`: Tries ENA first, falls back to SRA Toolkit if ENA links unavailable.

### Other fetch options

- `--out <SUBDIR>`: Output to specific subdirectory under `data/fastq/` (default: `data/fastq/`)
- `--parallel N`: Number of concurrent downloads (default: 4)
- `--threads N`: Number of SRA Toolkit threads (default: 4 or `THREADS` env var)

### After downloading FASTQs

FASTQs from GSE52778 do not come in correct naming format. Must use renaming helper script.
'''bash
bash scripts/rename_fastqs.sh
'''

Continue with the pipeline:

```bash
THREADS=8 bash salmon_pipeline.sh all
```

Then check Salmon mapping rates:

```bash
conda activate rnaseq
python scripts/check_salmon_qc.py --salmon_dir out/salmon --min_mapping 0.6
```

Finally, run DESeq2 analysis with the provided sample table:

```bash
conda activate rnaseq-r
Rscript scripts/run_deseq2.R \
  --quant_dir out/salmon \
  --gtf data/references/gtf/Homo_sapiens.GRCh38.115.gtf \
  --sample_table examples/sample_table.csv \
  --group_col condition \
  --project_name GSE52778_test \
  --padj_thresh 0.05 --lfc_thresh 0.5
```

## Expected Outputs

- After running the Salmon pipeline (`bash salmon_pipeline.sh all`):
  - `out/trimmed/`: Trimmed read pairs (`*_R1_trimmed.fastq.gz`, `*_R2_trimmed.fastq.gz`) from fastp.
  - `out/fastqc_trimmed/`: FastQC reports on trimmed reads (`*.html`, `*.zip`).
  - `out/multiqc/`: MultiQC summary (`multiqc_report.html`) aggregating fastp and trimmed FastQC reports.
  - `out/salmon/<sample>/`: Salmon quantification per sample (`quant.sf`, `lib_format_counts.json`, `meta_info.json`).
  - `logs/`: fastp JSON/HTML reports (`*.fastp.json`, `*.fastp.html`) with detailed trimming statistics.

- After running the QC check (`scripts/check_salmon_qc.py`):
  - Terminal output showing pass/fail status for each sample
  - Optional summary file if `--out` flag is used

- After running the DESeq2 wrapper (`scripts/run_deseq2.R`):
  - `out/deseq2/<PROJECT>_DESeq2_full_results.csv`: Per‑gene statistics (baseMean, log2FC, p-value, padj).
  - `out/deseq2/<PROJECT>_volcano_plot.pdf`: Volcano plot with thresholds.
  - `out/deseq2/<PROJECT>_PCA_plot.pdf`: PCA on VST-transformed counts.
  - `out/deseq2/<PROJECT>_raw_gene_counts.csv` and `_vst_norm_counts.csv`: Count matrices.
  - Enrichment outputs: Up/Down regulated enrichment Excel files. Note these are the enrichments on the top 200 genes within the specified l2FC & padj cutoffs.

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
  - `--project_name` prefix added to all outputs (e.g., `CU25_*.csv`)

### Rmd parameter: install_pkgs

- The Rmd has a parameter `install_pkgs` (default: false) that gates any `install.packages`/`BiocManager::install` calls.
- With the provided `environment-r.yml`, installs are not needed; leave `install_pkgs: false`.
- If running outside conda and you need the Rmd to install its own dependencies during render, set it in the YAML header or override at render time, for example:
  - Rscript -e "rmarkdown::render('tximport_deseq2.rmd', params=list(install_pkgs=TRUE))"

## Notes on practices and reproducibility

- This pipeline aligns with Bioconductor’s RNA-seq gene-level workflow guidance for design, dispersion, and multiple testing.
  - Bioconductor workflow: https://bioconductor.org/packages/release/workflows/html/rnaseqGene.html
- Please cite the tools used:
  - Salmon: Patro et al., Nature Methods 2017
  - DESeq2: Love et al., Genome Biology 2014
  - tximport: Soneson et al., F1000Research 2015
  - MultiQC: Ewels et al., Bioinformatics 2016
  - Ensembl/biomaRt for annotation

## Contributing

- Issues and PRs welcome for small, focused improvements (docs, minor fixes, portability). Please avoid changing the core analysis logic in `tximport_deseq2.rmd` unless requested.
- Keep changes minimal and backward compatible. For larger ideas, open an issue first to discuss scope.
- Style: keep bash scripts simple and echo clear progress; R changes should follow existing structure and use parameters where possible.

## Citations

- Salmon: Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods. 2017.
- fastp: Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018.
- DESeq2: Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014.
- tximport: Soneson C, Love MI, Robinson MD. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research. 2015.
- MultiQC: Ewels P et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016.
- Bioconductor workflow: "RNA-seq workflow: gene-level exploratory analysis and differential expression" (rnaseqGene).

## Repository layout

```
.
├── README.md, LICENSE, .gitignore
├── environment.yml (CLI dependencies)
├── environment-r.yml (R dependencies)
├── environment.lock.yml (pinned CLI versions)
├── environment-r.lock.yml (pinned R versions)
│
├── salmon_pipeline.sh (main pipeline orchestrator)
├── tximport_deseq2.rmd (DESeq2 analysis Rmarkdown)
│
├── scripts/
│   ├── get_refs.sh (download Ensembl references)
│   ├── fetch_test_fastqs.sh (fetch FASTQs from GEO/SRA)
│   ├── rename_fastqs.sh (standardize FASTQ naming)
│   ├── check_salmon_qc.py (verify Salmon mapping rates)
│   ├── run_deseq2.R (wrapper for DESeq2 analysis)
│   └── debug_deseq2.R (utility script)
│
├── data/
│   ├── fastq/ (place raw paired-end FASTQs here)
│   └── references/ (created by get_refs.sh)
│       ├── gtf/ (GTF files)
│       └── fa/ (FASTA cDNA files)
│
├── out/ (pipeline outputs)
│   ├── trimmed/
│   ├── fastqc_trimmed/
│   ├── multiqc/
│   ├── salmon/
│   └── deseq2/
│
├── logs/ (FastQC and fastp logs/reports)
│
└── examples/
    └── sample_table.csv (example metadata)
```

## Scripts and utilities

### scripts/fetch_test_fastqs.sh

Downloads FASTQs from NCBI GEO/SRA. Resolves GEO accessions (GSE) → SRA projects (SRP) → run accessions (SRR).

**Key features:**
- Automatic GEO accession resolution via NCBI GEO and ENA APIs
- Multiple download backends: SRA Toolkit, ENA (curl), or auto-fallback
- Parallel downloads with configurable thread/job count
- MD5 checksum verification (when available)
- Automatic gzip compression with `pigz` or `gzip`

**Usage examples:**
```bash
# From GEO accession (auto-resolves to SRR runs)
bash scripts/fetch_test_fastqs.sh --geo GSE52778

# Specific SRR runs
bash scripts/fetch_test_fastqs.sh --runs SRR1039508,SRR1039509,SRR1039512,SRR1039513

# Custom output directory
bash scripts/fetch_test_fastqs.sh --geo GSE52778 --out data/fastq/myrun

# Parallel downloads with ENA
bash scripts/fetch_test_fastqs.sh --geo GSE52778 --method ena --parallel 8
```

### scripts/get_refs.sh

Downloads Ensembl GTF and cDNA FASTA files for Salmon index building.

**Key features:**
- Supports human (GRCh38) and mouse (GRCm39) genomes
- Optional Ensembl release pinning (default: latest)
- Direct GTF flavor selection (default: auto-prefer plain GTF)
- Automatic decompression
- SHA256 checksums written to `data/references/README.md`

**Usage examples:**
```bash
# Human (default)
bash scripts/get_refs.sh --species human --build GRCh38

# Mouse
bash scripts/get_refs.sh --species mouse --build GRCm39

# Specific release
bash scripts/get_refs.sh --species human --build GRCh38 --release 110

# Explicit URLs
bash scripts/get_refs.sh --gtf_url https://... --fasta_url https://...
```

### scripts/rename_fastqs.sh

Standardizes FASTQ file naming from SRA-style (`SAMPLE_1.fastq.gz`, `SAMPLE_1.fq.gz`) to pipeline-expected format (`SAMPLE_R1_001.fastq.gz`).

**Key features:**
- Converts `_1/_2` to `_R1_001/_R2_001`
- Normalizes `.fq.gz` to `.fastq.gz`
- Dry-run mode to preview changes
- Custom directory support

**Usage examples:**
```bash
# Preview changes
bash scripts/rename_fastqs.sh --dry-run

# Apply changes
bash scripts/rename_fastqs.sh

# Custom folder
bash scripts/rename_fastqs.sh --dir /path/to/fastqs
```

### scripts/check_salmon_qc.py

Validates Salmon mapping rates against a minimum threshold. Useful for QC filtering before DESeq2.

**Key features:**
- Per-sample mapping rate reporting
- Pass/fail status based on threshold
- Optional summary file export
- Exit code 1 if any samples fail (suitable for CI/CD)

**Usage examples:**
```bash
# Check with 60% threshold
python scripts/check_salmon_qc.py --salmon_dir out/salmon --min_mapping 0.6

# Save summary
python scripts/check_salmon_qc.py --salmon_dir out/salmon --min_mapping 0.6 --out qc_summary.txt

# Strict threshold (80%)
python scripts/check_salmon_qc.py --salmon_dir out/salmon --min_mapping 0.8
```

### scripts/run_deseq2.R

Wrapper around `tximport_deseq2.rmd` for headless DESeq2 analysis. Handles GTF processing, tximport, and differential expression.

**Key features:**
- Automatic species detection (human vs mouse) for biomaRt annotation
- Gene symbol annotation from biomaRt (Ensembl → HGNC/MGI)
- All count matrices include gene symbols
- Volcano plots with labeled top genes
- Gene enrichment analysis (KEGG, GO, Reactome, etc.)
- Configurable p-value and log2FC thresholds

**Usage examples:**
```bash
# Basic human analysis
Rscript scripts/run_deseq2.R \
  --quant_dir out/salmon \
  --gtf data/references/gtf/Homo_sapiens.GRCh38.115.gtf \
  --sample_table examples/sample_table.csv \
  --group_col condition \
  --project_name MyProject

# Custom thresholds
Rscript scripts/run_deseq2.R \
  --quant_dir out/salmon \
  --gtf data/references/gtf/Homo_sapiens.GRCh38.115.gtf \
  --sample_table examples/sample_table.csv \
  --group_col condition \
  --project_name MyProject \
  --padj_thresh 0.01 \
  --lfc_thresh 1.0 \
  --out_dir results/deseq2
```

## Pipeline Design Rationale

- **fastp over Trimmomatic**: ~5x faster with automatic adapter detection and poly-G trimming
- **No raw FastQC**: fastp handles all QC issues; post-trim FastQC + Salmon metrics provide sufficient validation
- **MultiQC integration**: Automatically aggregates fastp, trimmed FastQC, and Salmon reports
- **Automatic gene annotation**: Species detection from GTF filename eliminates manual configuration
