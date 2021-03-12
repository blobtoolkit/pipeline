rule unchunk_blast:
    """
    Unchunk chunked blast results.
    """
    input:
        "{assembly}.diamond.reference_proteomes.out.raw"
    output:
        "{assembly}.diamond.reference_proteomes.out"
    params:
        max_target_seqs = config["similarity"]["max_target_seqs"]
    threads: 1
    log:
        "logs/{assembly}/unchunk_blast.log"
    benchmark:
        "logs/{assembly}/unchunk_blast.benchmark.txt"
    script:
        "../scripts/unchunk_blast.py"
