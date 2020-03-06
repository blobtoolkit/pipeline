rule blobtoolkit_add_cov:
    """
    Use BlobTools2 add to add coverage to a BlobDir from BAM files.
    """
    input:
        meta = "%s%s/identifiers.json" % (config['assembly']['prefix'], rev),
        bam = expand("%s.{sra}.bam" % asm, sra=list_sra_accessions(reads))
    output:
        expand("%s%s/{sra}_cov.json" % (config['assembly']['prefix'], rev), sra=list_sra_accessions(reads))
    params:
        id = "%s%s" % (config['assembly']['prefix'], rev),
        covs = lambda wc: ' --cov '.join(["%s.%s.bam=%s" % (config['assembly']['prefix'], sra, sra)
                                         for sra in list_sra_accessions(reads)])
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('blobtoolkit_add_cov', 1)
    log:
        "logs/%s/blobtoolkit_add_cov.log" % config['assembly']['prefix']
    benchmark:
        "logs/%s/blobtoolkit_add_cov.benchmark.txt" % config['assembly']['prefix']
    resources:
        threads = get_threads('blobtoolkit_add_cov', 1),
        btk = 1
    shell:
        'blobtools replace \
            --cov {params.covs} \
            --threads {threads} \
            {params.id} > {log} 2>&1'
