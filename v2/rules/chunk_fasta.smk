rule chunk_fasta:
    """
    Split long contigs into chunks.
    """
    input:
        fasta = "%s/{assembly}.windowmasker.fasta" % windowmasker_path,
    output:
        "{assembly}.fasta.chunks"
    params:
        chunk = set_blast_chunk(config),
        overlap = set_blast_chunk_overlap(config),
        max_chunks = set_blast_max_chunks(config),
        min_length = set_blast_min_length(config),
        fasta = lambda wc: "%s.fasta.bed" % wc.assembly
    threads: 1
    log:
        "logs/{assembly}/chunk_fasta.log"
    benchmark:
        "logs/{assembly}/chunk_fasta.benchmark.txt"
    script:
        "../scripts/chunk_fasta.py"
