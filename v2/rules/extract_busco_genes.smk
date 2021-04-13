rule extract_busco_genes:
    """
    Extract busco genes into a single fasta file.
    """
    input:
        busco = expand("%s/%s.busco.{lineage}/full_table.tsv.gz" % (busco_path, config["assembly"]["prefix"]), lineage=config['busco']['lineages']),
    output:
        "{assembly}.busco_genes.fasta"
    params:
        evalue = config["similarity"]["evalue"],
        max_target_seqs = config["similarity"]["max_target_seqs"],
        taxid = config["taxon"]["taxid"]
    threads: 32
    log:
        "logs/{assembly}/extract_busco_genes.log"
    benchmark:
        "logs/{assembly}/extract_busco_genes.benchmark.txt"
    shell:
        """(for file in {input.busco}; do \
            echo $file; \
        done) &> {log} && exit 1"""