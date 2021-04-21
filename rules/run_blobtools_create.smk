rule run_blobtools_create:
    """
    Run blobtools create.
    """
    input:
        length = "%s/%s.stats.length.bed" % (stats_path, config["assembly"]["prefix"]),
        gc = "%s/%s.stats.gc.bed" % (stats_path, config["assembly"]["prefix"]),
        busco = "%s/%s.busco.%s/full_table.tsv.gz" % (busco_path, config["assembly"]["prefix"], config['busco']['lineages'][0]),
        blastx = "%s/%s.diamond.reference_proteomes.out" % (diamond_path, config["assembly"]["prefix"]),
        blastn = "%s/%s.blastn.nt.out" % (blastn_path, config["assembly"]["prefix"]),
        taxdump = config["settings"]["taxdump"],
        yaml = "%s.meta.yaml" % config["assembly"]["prefix"]
    output:
        "%s/meta.json" % blobdir_name(config),
        "%s/identifiers.json" % blobdir_name(config),
        "%s/%s_phylum.json" % (blobdir_name(config), similarity_setting(config, "diamond_blastx", "taxrule")),
        expand("%s/{sra}_cov.json" % blobdir_name(config), sra=reads_by_prefix(config).keys()),
        "%s/%s_busco.json" % (blobdir_name(config), config['busco']['lineages'][0])
    params:
        statsdir = stats_path,
        busco = lambda wc: " --busco ".join(expand("%s/%s.busco.{lineage}/full_table.tsv.gz" % (busco_path, config["assembly"]["prefix"]), lineage=config['busco']['lineages'])),
        cov = blobtools_cov_flag(config),
        blobdir = blobdir_name(config),
        taxrule = similarity_setting(config, "diamond_blastx", "taxrule"),
        evalue = similarity_setting(config, "diamond_blastx", "import_evalue"),
        max_target_seqs = similarity_setting(config, "diamond_blastx", "import_max_target_seqs")
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
            --hits {input.blastx} \
            --hits {input.blastn} \
            --evalue {params.evalue} \
            --hit-count {params.max_target_seqs} \
            --threads {threads} \
            {params.blobdir} > {log} 2>&1"""
