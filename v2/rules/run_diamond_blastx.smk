rule run_diamond_blastx:
    """
    Run Diamond blastx to search protein database with assembly query.
    """
    input:
        fasta = "{assembly}.chunks.fasta",
        dmand = "reference_proteomes.dmnd"
    output:
        "{assembly}.diamond.reference_proteomes.out.raw"
    params:
        db = lambda wc: wc.diamond_db,
        evalue = lambda wc: similarity_config(config)["reference_proteomes"]["evalue"],
        max_target_seqs = lambda wc: similarity_config(config)["reference_proteomes"]["max_target_seqs"],
        taxid = lambda wc: config["taxonomy"]["taxid"]
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
            --taxon-exclude {params.taxon_id} \
            > {output} 2> {log}"""