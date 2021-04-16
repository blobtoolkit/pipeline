rule run_blobtools_create:
    """
    Run blobtools create.
    """
    input:
        length = "%s/%s.stats.length.bed" % (stats_path, config["assembly"]["prefix"]),
        gc = "%s/%s.stats.gc.bed" % (stats_path, config["assembly"]["prefix"]),
        busco = expand("%s/%s.busco.{lineage}/full_table.tsv.gz" % (busco_path, config["assembly"]["prefix"]), lineage=config['busco']['lineages']),
        blast = "%s/%s.diamond.reference_proteomes.out" % (diamond_path, config["assembly"]["prefix"]),
        taxdump = config["settings"]["taxdump"],
        yaml = "%s.meta.yaml" % config["assembly"]["prefix"]
    output:
        "%s/meta.json" % blobdir_name(config),
        "%s/identifiers.json" % blobdir_name(config),
        "%s/%s_phylum.json" % (blobdir_name(config), config["similarity"]["taxrule"]),
        expand("%s/{sra}_cov.json" % blobdir_name(config), sra=reads_by_prefix(config).keys()),
        expand("%s/{lineage}_busco.json" % blobdir_name(config), lineage=config['busco']['lineages'])
    params:
        statsdir = stats_path,
        busco = lambda wc: " --busco ".join(expand("%s/%s.busco.{lineage}/full_table.tsv.gz" % (busco_path, config["assembly"]["prefix"]), lineage=config['busco']['lineages'])),
        cov = blobtools_cov_flag(config),
        blobdir = blobdir_name(config),
        taxrule = config["similarity"]["taxrule"]
    threads: 4
    log:
        "logs/%s/run_blobtools_create.log" % config["assembly"]["prefix"]
    benchmark:
        "logs/%s/run_blobtools_create.benchmark.txt" % config["assembly"]["prefix"]
    shell:
        """blobtools replace \
            --beddir {params.statsdir} \
            --meta {input.yaml} \
            --taxdump {input.taxdump} \
            --taxrule {params.taxrule} \
            --busco {params.busco} \
            --hits {input.blast} \
            --threads {threads} \
            {params.blobdir} > {log} 2>&1"""
