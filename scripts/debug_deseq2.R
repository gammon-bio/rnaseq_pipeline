#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(readr)
  library(GenomicFeatures)
  library(txdbmaker)
  library(dplyr)
  library(tibble)
  library(AnnotationDbi)
})

# Parameters
quant_dir <- "out/salmon"
gtf_file <- "data/references/gtf/Homo_sapiens.GRCh38.115.gtf"
sample_table <- read.csv("examples/sample_table.csv", row.names = 1)

# Debug step 1: Check GTF parsing
message("Checking GTF parsing...")
txdb <- makeTxDbFromGFF(gtf_file)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
print("tx2gene structure:")
str(tx2gene)

# Debug step 2: Check Salmon files
message("\nChecking Salmon files...")
files <- file.path(quant_dir, row.names(sample_table), "quant.sf")
print("Salmon files:")
print(files)
print("Files exist:")
print(file.exists(files))

# Debug step 3: Import counts with version number handling
message("\nImporting counts...")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
print("txi structure:")
str(txi)

message("\nDone debugging.")
