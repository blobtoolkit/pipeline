rule chunk_fasta_by_busco:
    """
    Split long contigs into chunks containing busco genes.
    """
    input:
        fasta = "{assembly}.fasta",
        busco = "%s/{assembly}.busco.%s.tsv" % (busco_path, config["busco"]["lineages"][0])
    output:
        "{assembly}.fasta.chunks"
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
