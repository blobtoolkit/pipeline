rule fetch_assembly:
    """
    Fetch a remote assembly from EBI or NCBI.
    """
    output:
        fa = '{assembly}.fasta'
    params:
        url = lambda wc: prepare_ncbi_assembly_url(config['assembly']['accession'], config['assembly']['alias'])
    wildcard_constraints:
        assembly = r'\w+'
    conda:
        '../envs/fetch.yaml'
    threads: get_threads('fetch_assembly', 1)
    log:
        'logs/{assembly}/fetch_assembly.log'
    benchmark:
        'logs/{assembly}/fetch_assembly.benchmark.txt'
    resources:
        download = 1,
        threads = get_threads('fetch_assembly', 1)
    shell:
        '(curl -s {params.url} | \
        pigz -d > {output.fa}) 2> {log}'
