rule validate_dataset:
    """
    Run BlobToolKit validator on a dataset to check all expected fields are present.
    """
    input:
        cov = expand("%s%s/{sra}_cov.json" % (asm, rev), sra=list_sra_accessions(reads)),
        tax = "%s%s/%s_phylum_positions.json" % (asm, rev, taxrule_name()),
        busco = expand("%s%s/{lineage}_busco.json" % (asm, rev), lineage=config['busco']['lineages']),
        ids = "%s%s/identifiers.json" % (asm, rev)
    output:
        temp('{assembly}.valid')
    params:
        assembly = lambda wc: wc.assembly
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('validate_dataset', 1)
    log:
        'logs/{assembly}/validate_dataset.log'
    benchmark:
        'logs/{assembly}/validate_dataset.benchmark.txt'
    shell:
        'validate.py {params.assembly}/meta.json > {log} 2>&1 \
        && touch {params.assembly}.valid'
