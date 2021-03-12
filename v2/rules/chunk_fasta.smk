rule chunk_fasta:
    """
    Split long contigs into chunks.
    """
    input:
        fasta = "{assembly}.fasta"
    output:
        "{assembly}.fasta.chunks"
    params:
        chunk = 100000,
        overlap = 0,
        max_chunks = 10
    threads: 1
    log:
        "logs/{assembly}/chunk_fasta.log"
    benchmark:
        "logs/{assembly}/chunk_fasta.benchmark.txt"
    script:
        "../scripts/chunk_fasta.py"
