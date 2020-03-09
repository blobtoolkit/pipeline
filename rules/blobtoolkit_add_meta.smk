rule blobtoolkit_add_meta:
    """
    Use BlobTools2 replace to update metadata in a BlobDir.
    """
    input:
        yaml = '{assembly}.meta.yaml',
        lineages = "%s/taxidlineage.dmp" % config['settings']['taxonomy']
    output:
        temp('{assembly}.meta.updated')
    params:
        taxdump = taxdump_dir,
        id = lambda wc: "%s%s" % (wc.assembly, rev),
        taxid = config['taxon']['taxid'],
        path = config['settings']['blobtools2_path']
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('blobtoolkit_create', 1)
    log:
        'logs/{assembly}/blobtoolkit_create.log'
    benchmark:
        'logs/{assembly}/blobtoolkit_create.benchmark.txt'
    resources:
        threads = get_threads('blobtoolkit_create', 1),
        btk = 1
    shell:
        'PATH={params.path}:$PATH && \
        blobtools replace \
            --meta {input.yaml} \
            --taxid {params.taxid} \
            --taxdump {params.taxdump} \
            {params.id} > {log} 2>&1'
