#!/usr/bin/env Rscript

# run_deseq2.R — Wrapper to render tximport_deseq2.rmd with parameters
# Created: 2025-09-03
# Usage examples:
#   Rscript scripts/run_deseq2.R \
#     --quant_dir out/salmon \
#     --gtf data/references/gtf/Homo_sapiens.GRCh38.110.gtf \
#     --sample_table examples/sample_table.csv \
#     --group_col condition
#
#   Rscript scripts/run_deseq2.R \
#     --quant_dir out/salmon \
#     --tx2gene tx2gene.csv \
#     --sample_table examples/sample_table.csv
#
#   Rscript scripts/run_deseq2.R \
#     --tximport_rds path/to/txi.rds \
#     --sample_table examples/sample_table.csv

suppressPackageStartupMessages({
  library(optparse)
  library(rmarkdown)
})


option_list <- list(
  make_option(c("--quant_dir"), type = "character", default = NULL,
              help = "Directory with per-sample Salmon outputs (quant.sf)"),
  make_option(c("--gtf"), type = "character", default = NULL,
              help = "Path to GTF file (used by Rmd to build tx2gene)"),
  make_option(c("--tx2gene"), type = "character", default = NULL,
              help = "CSV with columns TXNAME,GENEID to use for tximport"),
  make_option(c("--tximport_rds"), type = "character", default = NULL,
              help = "Precomputed tximport object (.rds); skips re-import"),
  make_option(c("--sample_table"), type = "character", default = NULL,
              help = "CSV with sample metadata (rownames = sample IDs)"),
  make_option(c("--group_col"), type = "character", default = "condition",
              help = "Column in sample_table used for DE design [default %default]"),
  make_option(c("--padj_thresh"), type = "double", default = 0.05,
              help = "Adjusted P-value threshold [default %default]"),
  make_option(c("--lfc_thresh"), type = "double", default = 0.5,
              help = "Absolute log2FC threshold for volcano lines [default %default]"),
  make_option(c("--out_dir"), type = "character", default = "out/deseq2",
              help = "Output directory for DESeq2 artifacts [default %default]"),
  make_option(c("--project_name"), type = "character", default = "project",
              help = "Project name prefix for output files [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$sample_table)) {
  stop("--sample_table is required", call. = FALSE)
}

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# Prepare params for Rmd
params <- list(
  quant_dir   = if (!is.null(opt$quant_dir)) opt$quant_dir else "out/salmon",
  gtf_file    = if (!is.null(opt$gtf)) opt$gtf else "data/references/gtf/Mus_musculus.GRCm39.109.gtf",
  tximport_rds = NULL,
  sample_table = opt$sample_table,
  group_col    = opt$group_col,
  out_dir      = opt$out_dir,
  padj_thresh  = opt$padj_thresh,
  lfc_thresh   = opt$lfc_thresh,
  project_name = opt$project_name
)

# Option branch: tximport_rds provided directly
if (!is.null(opt$tximport_rds)) {
  params$tximport_rds <- opt$tximport_rds
} else if (!is.null(opt$quant_dir) && !is.null(opt$tx2gene)) {
  # Build tximport using provided tx2gene mapping, then pass RDS to Rmd
  message("[run_deseq2] Importing with tx2gene CSV: ", opt$tx2gene)
  suppressPackageStartupMessages({
    library(tximport)
    library(readr)
  })
  tx2gene <- readr::read_csv(opt$tx2gene, show_col_types = FALSE)
  if (!all(c("TXNAME","GENEID") %in% colnames(tx2gene))) {
    stop("--tx2gene must have columns: TXNAME, GENEID", call. = FALSE)
  }
  quant_dirs <- list.dirs(opt$quant_dir, recursive = FALSE, full.names = TRUE)
  files <- file.path(quant_dirs, "quant.sf")
  names(files) <- basename(quant_dirs)
  txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
  txi_path <- file.path(opt$out_dir, "txi.rds")
  saveRDS(txi, txi_path)
  params$tximport_rds <- txi_path
} else {
  # quant_dir + gtf path are passed through; Rmd will build tx2gene
  if (is.null(opt$quant_dir) || is.null(opt$gtf)) {
    message("[run_deseq2] Using Rmd defaults for quant_dir/gtf. You can pass --quant_dir and --gtf for explicit control.")
  }
}

# Render Rmd (resolve path relative to this script)
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub('^--file=', '', args_all[grep('^--file=', args_all)])
script_dir <- if (length(file_arg) == 1) dirname(normalizePath(file_arg)) else getwd()
rmd_path <- normalizePath(file.path(script_dir, "..", "tximport_deseq2.rmd"))

message("[run_deseq2] Rendering: ", rmd_path)
message("[run_deseq2] Params:")
param_strings <- sapply(names(params), function(n) {
  val <- if (is.null(params[[n]])) "NULL" else as.character(params[[n]])
  paste0("  ", n, " = ", val)
})
message(paste(param_strings, collapse = "\n"))

# Render with proper error handling and output visibility
tryCatch({
  rmarkdown::render(
    input = rmd_path,
    params = params,
    envir = new.env(),
    quiet = FALSE
  )
  message("[run_deseq2] Done. Outputs under: ", opt$out_dir)
}, error = function(e) {
  stop("[run_deseq2] ERROR during Rmd rendering: ", e$message, call. = FALSE)
})
