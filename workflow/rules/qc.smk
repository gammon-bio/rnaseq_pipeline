"""FastQC on trimmed reads, Salmon mapping-rate QC gate, and MultiQC summary
(mirrors salmon_pipeline.sh steps 5/7 + scripts/check_salmon_qc.py)."""


rule fastqc_trimmed:
    input:
        os.path.join(OUT, "trimmed", "{sample}_R{r}_trimmed.fastq.gz"),
    output:
        html=os.path.join(OUT, "fastqc_trimmed", "{sample}_R{r}_trimmed_fastqc.html"),
        zip=os.path.join(OUT, "fastqc_trimmed", "{sample}_R{r}_trimmed_fastqc.zip"),
    threads: 1
    conda:
        "../envs/fastqc.yaml"
    log:
        os.path.join("logs", "fastqc_{sample}_R{r}.log"),
    shell:
        "fastqc -t {threads} -o $(dirname {output.html}) {input} 2>&1 | tee {log}"


rule salmon_qc_check:
    input:
        quant_files(),
    output:
        QC_SUMMARY,
    params:
        salmon_dir=os.path.join(OUT, "salmon"),
        min_mapping=config["min_mapping_rate"],
    conda:
        "../envs/multiqc.yaml"
    log:
        os.path.join("logs", "salmon_qc_check.log"),
    shell:
        "python scripts/check_salmon_qc.py --salmon_dir {params.salmon_dir} "
        "--min_mapping {params.min_mapping} --out {output} 2>&1 | tee {log}"


rule multiqc:
    input:
        fastqc=fastqc_files(),
        quant=quant_files(),
    output:
        MULTIQC_REPORT,
    params:
        scan=[os.path.join(OUT, "fastqc_trimmed"), os.path.join(OUT, "salmon"), "logs"],
        outdir=os.path.join(OUT, "multiqc"),
    conda:
        "../envs/multiqc.yaml"
    log:
        os.path.join("logs", "multiqc.log"),
    shell:
        "multiqc {params.scan} -o {params.outdir} -n multiqc_report.html "
        "--force 2>&1 | tee {log}"
