rule unchunk_blastn:
    """
    Unchunk chunked blastn results.
    """
    input:
        "{assembly}.blastn.nt.out.raw"
    output:
        "{assembly}.blastn.nt.out"
    params:
        max_target_seqs = similarity_setting(config, "blastn", "max_target_seqs")
    threads: 1
    log:
        "logs/{assembly}/unchunk_blastn.log"
    benchmark:
        "logs/{assembly}/unchunk_blastn.benchmark.txt"
    script:
        "../lib/unchunk_blast.py"
