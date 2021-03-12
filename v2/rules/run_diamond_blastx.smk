rule run_diamond_blastx:
    """
    Run Diamond blastx to search protein database with assembly query.
    """
    input:
        fasta = "{assembly}.fasta.chunks",
        dmnd = "%s/%s.dmnd" % (config["similarity"]["path"], config["similarity"]["name"])
    output:
        "{assembly}.diamond.reference_proteomes.out.raw"
    params:
        db = config["similarity"]["path"],
        evalue = config["similarity"]["evalue"],
        max_target_seqs = config["similarity"]["max_target_seqs"],
        taxid = config["taxonomy"]["taxid"]
    threads: 32
    log:
        "logs/{assembly}/run_diamond_blastx.log"
    benchmark:
        "logs/{assembly}/run_diamond_blastx.benchmark.txt"
    shell:
        """diamond blastx \
            --query {input.fasta} \
            --db {params.db} \
            --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --max-target-seqs {params.max_target_seqs} \
            --max-hsps 1 \
            --evalue {params.evalue} \
            --threads {threads} \
            --taxon-exclude {params.taxid} \
            > {output} 2> {log}"""