rule run_blobtools_add:
    """
    Run blobtools add.
    """
    input:
        meta = "%s/meta.json" % blobdir_name(config),
        blastp = "%s/%s.diamond.busco_genes.out" % (diamond_blastp_path, config["assembly"]["prefix"]),
        taxdump = config["settings"]["taxdump"],
    output:
        "%s/busco_phylum.json" % blobdir_name(config)
    params:
        blobdir = blobdir_name(config),
        evalue = similarity_setting(config, "diamond_blastp", "import_evalue"),
        max_target_seqs = similarity_setting(config, "diamond_blastp", "import_max_target_seqs")
    threads: 4
    log:
        "logs/%s/run_blobtools_add.log" % config["assembly"]["prefix"]
    benchmark:
        "logs/%s/run_blobtools_add.benchmark.txt" % config["assembly"]["prefix"]
    shell:
        """blobtools replace \
            --taxdump {input.taxdump} \
            --taxrule blastp=busco \
            --hits {input.blastp} \
            --evalue {params.evalue} \
            --hit-count {params.max_target_seqs} \
            --threads {threads} \
            {params.blobdir} > {log} 2>&1"""
