rule chunk_nohit_fasta:
    """
    Split long contigs into chunks.
    """
    input:
        fasta = "{assembly}.nohit.fasta",
    output:
        fasta = "{assembly}.nohit.fasta.chunks"
    params:
        chunk = set_blast_chunk(config),
        overlap = set_blast_chunk_overlap(config),
        max_chunks = set_blast_max_chunks(config),
        min_length = set_blast_min_length(config)
    threads: 1
    log:
        "logs/{assembly}/chunk_fasta.log"
    benchmark:
        "logs/{assembly}/chunk_fasta.benchmark.txt"
    script:
        "../lib/chunk_fasta.py"
