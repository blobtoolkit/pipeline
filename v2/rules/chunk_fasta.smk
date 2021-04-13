rule chunk_fasta:
    """
    Split long contigs into chunks.
    """
    input:
        fasta = "%s/{assembly}.windowmasker.fasta" % windowmasker_path,
    output:
        "{assembly}.fasta.chunks"
    params:
        chunk = 100000,
        overlap = 0,
        max_chunks = 10,
        min_length = 1000
    threads: 1
    log:
        "logs/{assembly}/chunk_fasta.log"
    benchmark:
        "logs/{assembly}/chunk_fasta.benchmark.txt"
    script:
        "../scripts/chunk_fasta.py"
