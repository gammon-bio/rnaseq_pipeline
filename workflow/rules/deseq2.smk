"""tximport + DESeq2 differential expression.

Wraps the existing scripts/run_deseq2.R (-> tximport_deseq2.rmd) unchanged, so the
full multi-factor / LRT feature set is preserved. All DE settings come from the
config `deseq2:` block; optional flags (design/test/reduced/ref_level/contrast) are
assembled by deseq2_optional_args() and omitted when null.
"""


rule deseq2:
    input:
        quant=quant_files(),
        sample_table=config["deseq2"]["sample_table"],
        gtf=config["gtf"],
    output:
        deseq2_target(),
    params:
        quant_dir=os.path.join(OUT, "salmon"),
        out_dir=os.path.join(OUT, "deseq2"),
        group_col=config["deseq2"]["group_col"],
        project_name=config["deseq2"]["project_name"],
        padj=config["deseq2"]["padj_thresh"],
        lfc=config["deseq2"]["lfc_thresh"],
        optional=deseq2_optional_args,
    conda:
        "../../environment-r.yml"
    log:
        os.path.join("logs", "deseq2.log"),
    shell:
        "Rscript scripts/run_deseq2.R "
        "--quant_dir {params.quant_dir} --gtf {input.gtf} "
        "--sample_table {input.sample_table} --group_col {params.group_col} "
        "--project_name {params.project_name} --out_dir {params.out_dir} "
        "--padj_thresh {params.padj} --lfc_thresh {params.lfc} "
        "{params.optional} 2>&1 | tee {log}"
