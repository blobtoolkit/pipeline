rule blobtoolkit_remove_busco:
    """
    Remove BUSCO results from a BlobDir.
    """
    input:
        "{assembly}%s/meta.json" % rev
    output:
        temp("{assembly}%s.busco.removed" % rev)
    params:
        id = "{assembly}%s" % rev,
        path = config['settings']['blobtools2_path']
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('blobtoolkit_remove_busco', 1)
    log:
        'logs/{assembly}/blobtoolkit_remove_busco.log'
    benchmark:
        'logs/{assembly}/blobtoolkit_remove_busco.benchmark.txt'
    resources:
        btk = 1
    priority:
        200
    shell:
        'PATH={params.path}:$PATH && \
        blobtools remove \
            --busco \
            {params.id} > {log} 2>&1'
