rule blobtoolkit_remove_hits:
    """
    Remove similarity search results from a BlobDir.
    """
    input:
        "{assembly}%s/meta.json" % rev
    output:
        touch(temp("{assembly}%s.hits.removed" % rev))
    params:
        id = "{assembly}%s" % rev,
        path = config['settings']['blobtools2_path']
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('blobtoolkit_remove_hits', 1)
    log:
        'logs/{assembly}/blobtoolkit_remove_hits.log'
    benchmark:
        'logs/{assembly}/blobtoolkit_remove_hits.benchmark.txt'
    resources:
        btk = 1
    priority:
        200
    shell:
        'PATH={params.path}:$PATH && \
        blobtools remove \
            --hits \
            {params.id} > {log} 2>&1'
