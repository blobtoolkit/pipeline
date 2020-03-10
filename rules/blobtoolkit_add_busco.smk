rule blobtoolkit_add_busco:
    """
    import BUSCO results into BlobDir.
    """
    input:
        meta = "%s%s/identifiers.json" % (config['assembly']['prefix'], rev),
        tsv = expand("%s.busco.{lineage}.tsv" % config['assembly']['prefix'], lineage=config['busco']['lineages'])
    output:
        expand("%s%s/{lineage}_busco.json" % (config['assembly']['prefix'], rev), lineage=config['busco']['lineages'])
    params:
        id = "%s%s" % (config['assembly']['prefix'], rev),
        busco = ' --busco '.join(["%s.busco.%s.tsv" % (config['assembly']['prefix'], lineage)
                                  for lineage in config['busco']['lineages']]),
        path = config['settings']['blobtools2_path']
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('blobtoolkit_add_busco', 1)
    log:
        "logs/%s/blobtoolkit_add_busco.log" % config['assembly']['prefix']
    benchmark:
        "logs/%s/blobtoolkit_add_busco.benchmark.txt" % config['assembly']['prefix']
    resources:
        btk = 1
    shell:
        'PATH={params.path}:$PATH && \
        blobtools replace \
            --busco {params.busco} \
            {params.id} > {log} 2>&1'
