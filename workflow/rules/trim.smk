"""fastp adapter/quality trimming (mirrors salmon_pipeline.sh step 5)."""


rule fastp:
    input:
        r1=os.path.join(FASTQ_DIR, "{sample}_R1_001.fastq.gz"),
        r2=os.path.join(FASTQ_DIR, "{sample}_R2_001.fastq.gz"),
    output:
        r1=os.path.join(OUT, "trimmed", "{sample}_R1_trimmed.fastq.gz"),
        r2=os.path.join(OUT, "trimmed", "{sample}_R2_trimmed.fastq.gz"),
        json=os.path.join("logs", "{sample}.fastp.json"),
        html=os.path.join("logs", "{sample}.fastp.html"),
    params:
        extra=config["fastp_extra"],
    threads: THREADS
    conda:
        "../envs/fastp.yaml"
    log:
        os.path.join("logs", "fastp_{sample}.log"),
    shell:
        "fastp -i {input.r1} -I {input.r2} "
        "-o {output.r1} -O {output.r2} "
        "{params.extra} --thread {threads} "
        "--json {output.json} --html {output.html} 2>&1 | tee {log}"
