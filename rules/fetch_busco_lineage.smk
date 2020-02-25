rule fetch_busco_lineage:
    """
    Fetch BUSCO lineages
    """
    output:
        gz = config['busco']['lineage_dir']+'/{lineage}.tar.gz',
        cfg = config['busco']['lineage_dir']+'/{lineage}/dataset.cfg'
    params:
        lineage = lambda wc: wc.lineage,
        dir = lambda wc: config['busco']['lineage_dir']
    conda:
        '../envs/fetch.yaml'
    threads: get_threads('fetch_busco_lineage', 1)
    log:
        'logs/fetch_busco_lineage/{lineage}.log'
    benchmark:
        'logs/fetch_busco_lineage/{lineage}.benchmark.txt'
    resources:
        download = 1,
        threads = get_threads('fetch_busco_lineage', 1)
    shell:
        '(wget -q -O {output.gz} "https://busco.ezlab.org/datasets/{params.lineage}.tar.gz" \
        && tar xf {output.gz} -C {params.dir}) 2> {log}'
