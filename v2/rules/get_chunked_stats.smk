rule get_chunked_stats:
    """
    Get chunked sequence stats.
    """
    input:
        fasta = "%s/{assembly}.windowmasker.fasta" % windowmasker_path,
    output:
        mask = "{assembly}.stats.mask.bed",
        gc = "{assembly}.stats.gc.bed",
        gc_windows = "{assembly}.stats.gc_windows.bed",
        n = "{assembly}.stats.n.bed",
        n_windows = "{assembly}.stats.n_windows.bed",
        masked = "{assembly}.stats.masked.bed",
        masked_windows = "{assembly}.stats.masked_windows.bed",
        length = "{assembly}.stats.length.bed"
    params:
        chunk = set_blast_chunk(config),
        overlap = set_blast_chunk_overlap(config),
        max_chunks = set_blast_max_chunks(config),
        min_length = set_blast_min_length(config),
        bed = lambda wc: "%s.stats.bed" % wc.assembly
    threads: 1
    log:
        "logs/{assembly}/get_seq_stats.log"
    benchmark:
        "logs/{assembly}/get_seq_stats.benchmark.txt"
    script:
        "../scripts/chunk_fasta.py"
