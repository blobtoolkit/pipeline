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
        """(> {output}; \
        for TABLE in {input.busco}; do \
            if [ -s $TABLE ]; then \
                SEQS=${{TABLE/full_table.tsv.gz/busco_sequences.tar.gz}};
                tar xf $SEQS \
                    --to-command='FILE=$(basename $TAR_FILENAME); \
                                  awk -v busco=${{FILE%.faa}} '"'"'{{if($1 ~ /^>/){{print $1 "=" busco}} else {{print $1}}}}'"'"''; \
            fi; \
        done) >> {output} 2> {log}"""
