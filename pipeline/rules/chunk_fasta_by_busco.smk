rule chunk_fasta_by_busco:
    """
    Split long contigs into chunks containing busco genes.
    """
    input:
        fasta = "%s/{assembly}.windowmasker.fasta" % windowmasker_path,
        busco = "%s/{assembly}.busco.%s/full_table.tsv.gz" % (busco_path, config["busco"]["lineages"][0])
    output:
        fasta = "{assembly}.fasta.chunks"
    params:
        chunk = set_blast_chunk(config),
        overlap = set_blast_chunk_overlap(config),
        max_chunks = set_blast_max_chunks(config),
        min_length = set_blast_min_length(config)
    threads: 1
    log:
        "logs/{assembly}/chunk_fasta_by_busco.log"
    benchmark:
        "logs/{assembly}/chunk_fasta_by_busco.benchmark.txt"
    script:
        "../lib/chunk_fasta.py"
