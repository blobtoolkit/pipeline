rule generate_summary:
    """
    Use BlobTools2 filter to generate a dataset summary.
    """
    input:
        '{assembly}/cumulative.png'
    output:
        '{assembly}/summary.json'
    params:
        assembly = lambda wc: wc.assembly
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('generate_summary', 1)
    log:
        'logs/{assembly}/generate_summary.log'
    benchmark:
        'logs/{assembly}/generate_summary.benchmark.txt'
    shell:
        'blobtools filter --summary {output} {params.assembly} 2> {log}'
