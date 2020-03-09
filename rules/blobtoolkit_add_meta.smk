rule blobtoolkit_add_meta:
    """
    Use BlobTools2 replace to update metadata in a BlobDir.
    """
    input:
        yaml = '{assembly}.meta.yaml'
    output:
        temp('{assembly}.meta.updated')
    params:
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
            {params.id} > {log} 2>&1'
