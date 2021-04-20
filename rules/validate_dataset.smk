rule validate_dataset:
    """
    Run BlobToolKit validator on a dataset to check all expected fields are present.
    """
    input:
        copied = "{blobdir}.copied",
        cov = expand("{{blobdir}}/{sra}_cov.json", sra=reads_by_prefix(config).keys()),
        tax = "{blobdir}/%s_phylum_positions.json" % similarity_setting(config, "diamond_blastx", "taxrule"),
        busco = expand("{{blobdir}}/{lineage}_busco.json", lineage=config["busco"]["lineages"]),
        ids = "{blobdir}/identifiers.json"
    output:
        touch(temp("{blobdir}.valid"))
    params:
        blobdir = lambda wc: wc.blobdir
    threads: 1
    log:
        "logs/{blobdir}/validate_dataset.log"
    benchmark:
        "logs/{blobdir}/validate_dataset.benchmark.txt"
    shell:
        """validate.py {params.blobdir}/meta.json > {log} 2>&1"""
