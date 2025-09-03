# Data directory

This folder holds large, user-provided inputs and downloaded references. It is intentionally git-ignored (except for this README and .gitkeep files).

- fastq/
  - Place raw paired-end FASTQs here.
  - Expected names: <sample>_R1_001.fastq.gz and <sample>_R2_001.fastq.gz.
  - The pipeline reads from data/fastq/ and writes outputs under out/.

- references/
  - Populated by bash scripts/get_refs.sh.
  - Subfolders: gtf/ (annotation), fa/ (cDNA FASTA).
  - Script records exact URLs and SHA256 checksums in data/references/README.md.

Notes
- Do not commit large sequence files to git.
- If you maintain multiple genomes/releases, keep them in separate subfolders or rename files accordingly.
