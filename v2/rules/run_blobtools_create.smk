rule run_blobtools_create:
    """
    Run blobtools create.
    """
    input:
        fasta = "%s.fasta" % config["assembly"]["prefix"],
        busco = expand("%s.busco.{lineage}.tsv" % config["assembly"]["prefix"], lineage=config['busco']['lineages']),
        bam = expand("%s.{sra}.bam" % config["assembly"]["prefix"], sra=reads_by_prefix(config).keys()),
        blast = "%s.%s.out" % (config["assembly"]["prefix"], diamond_db_name(config)),
        taxdump = config["similarity"],
        yaml = "%s.meta.yaml" % config["assembly"]["prefix"]
    output:
        "%s/meta.json" % blobdir_name(config),
        "%s/identifiers.json" % blobdir_name(config),
        "%s/%s_phylum.json" % (blobdir_name(config), config["similarity"]["taxrule"]),
        expand("%s/{sra}_cov.json" % blobdir_name(config), sra=reads_by_prefix(config).keys()),
        expand("%s/{lineage}_busco.json" % blobdir_name(config), lineage=config['busco']['lineages'])
    params:
        busco = lambda wc: " --busco ".join(input.busco),
        bam = lambda wc: " --cov ".join(input.bam),
    threads: 30
    log:
        "logs/%s/run_blobtools_create.log" % config["assembly"]["prefix"]
    benchmark:
        "logs/%s/run_blobtools_create.benchmark.txt" % config["assembly"]["prefix"]
    shell:
        """blobtools replace \
            --fasta {input.fasta} \
            --meta {input.yaml} \
            --taxdump {input.taxdump} \
            --busco {params.busco} \
            --cov {params.cov} \
            --hits {input.blast}
            {params.id} > {log} 2>&1"""