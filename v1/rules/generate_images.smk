rule generate_images:
    """
    Use BlobTools2 view to generate a set of static images.
    """
    input:
        valid = '{assembly}.valid',
        cov = expand("%s%s/{sra}_cov.json" % (asm, rev), sra=list_sra_accessions(reads))
    output:
        '{assembly}/cumulative.png'
    params:
        assembly = lambda wc: wc.assembly,
        host = 'http://localhost',
        ports = '8000-8099'
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('generate_images', 3)
    log:
        'logs/{assembly}/generate_images.log'
    benchmark:
        'logs/{assembly}/generate_images.benchmark.txt'
    script:
        '../scripts/generate_static_images.py'
