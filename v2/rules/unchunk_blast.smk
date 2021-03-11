rule unchunk_blast:
    """
    Unchunk chunked blast results.
    """
    input:
        "{assembly}.diamond.reference_proteomes.out.raw"
    output:
        "{assembly}.diamond.reference_proteomes.out"
    params:
        max_target_seqs = lambda wc: similarity_config(config)["reference_proteomes"]["max_target_seqs"]
    threads: 1
    log:
        "logs/{assembly}/unchunk_blast.log"
    benchmark:
        "logs/{assembly}/unchunk_blast.benchmark.txt"
    script:
        "../scripts/unchunk_blast.py"
