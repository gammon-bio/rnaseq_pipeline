"""Optional reference download (wraps scripts/get_refs.sh).

DECOUPLED from the default DAG on purpose: the processing rules consume the
reference files (gentrome/decoys/transcriptome/gtf) as pre-existing inputs, so a
plain `snakemake all` never re-downloads references or rebuilds the protected
index. Run this explicitly when you need to (re)fetch a reference set:

    snakemake get_refs                                  # uses config refs: block
    snakemake get_refs --configfile config/config.test.yaml   # human, cDNA-only

The real reference files land under data/references/ as a side effect; this rule
only writes a marker so it stays out of the main dependency graph.
"""

localrules: get_refs


rule get_refs:
    output:
        touch(os.path.join("resources", "markers", "get_refs.done")),
    params:
        species=config["refs"]["species"],
        build=config["refs"]["build"],
        release=config["refs"]["release"],
        decoy="--decoy" if config["refs"].get("decoy", True) else "--no-decoy",
    log:
        os.path.join("logs", "get_refs.log"),
    shell:
        "bash scripts/get_refs.sh --species {params.species} --build {params.build} "
        "--release {params.release} {params.decoy} 2>&1 | tee {log}"
