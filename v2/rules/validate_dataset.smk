rule validate_dataset:
    """
    Run BlobToolKit validator on a dataset to check all expected fields are present.
    """
    input:
        cov = expand("{{blobdir}}/{sra}_cov.json", sra=list_sra_accessions(reads)),
        tax = expand("{{blobdir}}/{taxrule}_phylum_positions.json" % (asm, rev), taxrule=taxrule_name()),
        busco = expand("{{blobdir}}/{lineage}_busco.json" % (asm, rev), lineage=config["busco"]["lineages"]),
        ids = "{{blobdir}}/identifiers.json" % (asm, rev)
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
