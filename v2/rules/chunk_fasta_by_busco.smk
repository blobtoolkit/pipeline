rule chunk_fasta_by_busco:
    """
    Split long contigs into chunks containing busco genes.
    """
    input:
        fasta = "%s/{assembly}.windowmasker.fasta" % windowmasker_path,
        busco = "%s/{assembly}.busco.%s/full_table.tsv.gz" % (busco_path, config["busco"]["lineages"][0])
    output:
        "{assembly}.fasta.chunks"
    params:
        chunk = 100000,
        overlap = 0,
        max_chunks = 10,
        min_length = 1000
    threads: 1
    log:
        "logs/{assembly}/chunk_fasta_by_busco.log"
    benchmark:
        "logs/{assembly}/chunk_fasta_by_busco.benchmark.txt"
    script:
        "../scripts/chunk_fasta.py"
