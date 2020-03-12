rule blobtoolkit_remove_cov:
    """
    Remove read coverage from a BlobDir.
    """
    input:
        "{assembly}%s/meta.json" % rev
    output:
        temp("{assembly}%s.coverage.removed" % rev)
    params:
        id = "{assembly}%s" % rev,
        path = config['settings']['blobtools2_path']
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('blobtoolkit_remove_cov', 1)
    log:
        'logs/{assembly}/blobtoolkit_remove_cov.log'
    benchmark:
        'logs/{assembly}/blobtoolkit_remove_cov.benchmark.txt'
    resources:
        btk = 1
    priority:
        200
    shell:
        'PATH={params.path}:$PATH && \
        blobtools remove \
            --cov \
            {params.id} > {log} 2>&1'
