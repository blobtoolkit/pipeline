rule fetch_busco_lineage:
    """
    Fetch BUSCO lineages
    """
    output:
        gz = config['busco']['lineage_dir']+'/{lineage}.tar.gz',
        cfg = config['busco']['lineage_dir']+'/{lineage}/dataset.cfg'
    params:
        lineage = lambda wc: wc.lineage,
        dir = busco_dir,
        path = lambda wc: 'busco-archive.ezlab.org/v3/datasets' if wc.lineage.endswith('9') else 'busco-data.ezlab.org/v4/data/lineages',
        date = lambda wc: '' if wc.lineage.endswith('9') else '.2019-11-20'
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
        '(wget -q -O {params.dir}/{params.lineage}.tar.gz "https://{params.path}/{params.lineage}{params.date}.tar.gz" \
        && tar xf {params.dir}/{params.lineage}.tar.gz -C {params.dir}) 2> {log}'
