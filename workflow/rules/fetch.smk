"""Optional test-FASTQ download (wraps scripts/fetch_test_fastqs.sh + rename_fastqs.sh).

DECOUPLED from the default DAG (same rationale as refs.smk). Run explicitly to
populate data/fastq/ for the GSE52778 smoke test, then run the pipeline:

    snakemake fetch_test_data --configfile config/config.test.yaml
    snakemake all            --configfile config/config.test.yaml --sdm conda --conda-frontend mamba -c 8

It downloads the configured runs and then standardizes the filenames to
<sample>_R1_001.fastq.gz / _R2_001.fastq.gz so sample discovery picks them up on
the next invocation (glob_wildcards is evaluated when the workflow loads).
"""

localrules: fetch_test_data


rule fetch_test_data:
    output:
        touch(os.path.join("resources", "markers", "fetch_test_data.done")),
    params:
        geo=config["fetch"]["geo"],
        runs=config["fetch"]["runs"],
        method=config["fetch"]["method"],
    log:
        os.path.join("logs", "fetch_test_data.log"),
    shell:
        "bash scripts/fetch_test_fastqs.sh --geo {params.geo} --runs {params.runs} "
        "--method {params.method} 2>&1 | tee {log} && "
        "bash scripts/rename_fastqs.sh 2>&1 | tee -a {log}"
