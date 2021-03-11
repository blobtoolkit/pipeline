import os

rule add_summary_to_metadata:
    input:
        stats = expand("{{assembly}}.{sra}.bam.stats", sra=reads_by_prefix(config).keys())
    output:
        "{assembly}.meta.yaml"
    params:
        gitdir = git_dir
    threads: 1
    log:
        "logs/{assembly}/generate_metadata.log"
    benchmark:
        "logs/{assembly}/generate_metadata.benchmark.txt"
    script:
        "../scripts/generate_metadata.py"
