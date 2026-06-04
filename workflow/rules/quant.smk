"""Salmon index + quant (mirrors salmon_pipeline.sh step 6).

The salmon_index rule is built to NEVER rebuild an existing index:
  * inputs are ancient()  -> a refreshed reference timestamp can't trigger a rebuild
  * output is protected()  -> Snakemake refuses to delete/overwrite it
Register the prebuilt mouse index once with `snakemake --touch salmon_index`; thereafter
a normal run sees the directory present and skips this rule entirely.
"""


rule salmon_index:
    input:
        unpack(salmon_index_inputs),
    output:
        protected(directory(SALMON_INDEX)),
    params:
        decoy_flag=lambda wildcards, input: (
            "-d {}".format(input.decoys) if config.get("decoy_aware", True) else ""
        ),
    threads: THREADS
    conda:
        "../envs/salmon.yaml"
    log:
        os.path.join("logs", "salmon_index.log"),
    shell:
        "salmon index -t {input.ref} {params.decoy_flag} "
        "-i {output} -k 31 -p {threads} 2>&1 | tee {log}"


rule salmon_quant:
    input:
        r1=os.path.join(OUT, "trimmed", "{sample}_R1_trimmed.fastq.gz"),
        r2=os.path.join(OUT, "trimmed", "{sample}_R2_trimmed.fastq.gz"),
        index=ancient(SALMON_INDEX),
    output:
        quant=os.path.join(OUT, "salmon", "{sample}", "quant.sf"),
        meta=os.path.join(OUT, "salmon", "{sample}", "aux_info", "meta_info.json"),
    params:
        libtype=config["salmon_libtype"],
        extra=config["salmon_extra"],
    threads: THREADS
    conda:
        "../envs/salmon.yaml"
    log:
        os.path.join("logs", "salmon_quant_{sample}.log"),
    shell:
        "salmon quant -i {input.index} -l {params.libtype} "
        "-1 {input.r1} -2 {input.r2} -p {threads} {params.extra} "
        "-o $(dirname {output.quant}) 2>&1 | tee {log}"
