rule generate_summary:
    """
    Use BlobTools2 filter to generate a dataset summary.
    """
    input:
        '{assembly}/cumulative.png'
    output:
        '{assembly}/summary.json'
    params:
        assembly = lambda wc: wc.assembly,
        taxrule = taxrule_name()[0]
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('generate_summary', 1)
    log:
        'logs/{assembly}/generate_summary.log'
    benchmark:
        'logs/{assembly}/generate_summary.benchmark.txt'
    shell:
        'blobtools filter --summary {output} --taxrule {params.taxrule} {params.assembly} 2> {log}'
