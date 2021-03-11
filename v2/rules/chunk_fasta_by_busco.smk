rule chunk_fasta_by_busco:
    """
    Split long contigs into chunks containing busco genes.
    """
    input:
        fasta = "{assembly}.fasta",
        busco = "%s/{assembly}.busco.%s.tsv" % (busco_dir, config["busco"]["lineages"])
    output:
        "{assembly}.chunks.fasta"
    params:
        chunk = 100000,
        overlap = 500,
        max_chunks = 10
    threads: 1
    log:
        "logs/{assembly}/chunk_fasta_by_busco.log"
    benchmark:
        "logs/{assembly}/chunk_fasta_by_busco.benchmark.txt"
    script:
        "../scripts/chunk_fasta.py"
