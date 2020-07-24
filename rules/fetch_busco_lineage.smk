import re

def busco_date(wc):
    """find the date for a BUSCO lineage."""
    if wc.lineage.endswith('9'):
        return ''
    with open('busco_downloads/file_versions.tsv') as fh:
        for line in fh:
            if line.startswith(wc.lineage):
                parts = re.split(r'\t', line)
                return '.'+parts[1]
    return ''


rule fetch_busco_versions:
    """
    Fetch BUSCO versions
    """
    output:
        'busco_downloads/file_versions.tsv'
    conda:
        '../envs/fetch.yaml'
    threads: get_threads('fetch_busco_versions', 1)
    log:
        'logs/fetch_busco_versions.log'
    benchmark:
        'logs/fetch_busco_versions.benchmark.txt'
    resources:
        download = 1
    shell:
        '(wget -q -O {output} "https://busco-data.ezlab.org/v4/data/file_versions.tsv") 2> {log}'


rule fetch_busco_lineage:
    """
    Fetch BUSCO lineages
    """
    input:
        'busco_downloads/file_versions.tsv'
    output:
        gz = config['busco']['lineage_dir']+'/{lineage}.tar.gz',
        cfg = config['busco']['lineage_dir']+'/{lineage}/dataset.cfg'
    params:
        lineage = lambda wc: wc.lineage,
        dir = busco_dir,
        path = lambda wc: 'busco-archive.ezlab.org/v3/datasets' if wc.lineage.endswith('9') else 'busco-data.ezlab.org/v4/data/lineages',
        date = busco_date
    conda:
        '../envs/fetch.yaml'
    threads: get_threads('fetch_busco_lineage', 1)
    log:
        'logs/fetch_busco_lineage/{lineage}.log'
    benchmark:
        'logs/fetch_busco_lineage/{lineage}.benchmark.txt'
    resources:
        download = 1
    shell:
        '(wget -q -O {params.dir}/{params.lineage}.tar.gz "https://{params.path}/{params.lineage}{params.date}.tar.gz" \
        && tar xf {params.dir}/{params.lineage}.tar.gz -C {params.dir}) 2> {log}'
