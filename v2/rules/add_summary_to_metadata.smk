rule add_summary_to_metadata:
    input:
        stats = expand("{{assembly}}.{sra}.bam.stats", sra=reads_by_prefix(config).keys())
    output:
        "{assembly}.meta.yaml"
    threads: 1
    log:
        "logs/{assembly}/add_summary_to_metadata.log"
    benchmark:
        "logs/{assembly}/add_summary_to_metadata.benchmark.txt"
    script:
        "../scripts/add_summary_to_metadata.py"
