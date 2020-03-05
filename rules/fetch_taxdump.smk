rule fetch_taxdump:
    """
    Fetch NCBI taxonomy "names.dmp" and "nodes.dmp" files
    """
    output:
        "%s/names.dmp" % config['settings']['taxonomy'],
        "%s/nodes.dmp" % config['settings']['taxonomy'],
        "%s/taxidlineage.dmp" % config['settings']['taxonomy']
    params:
        outdir=taxdump_dir
    conda:
        '../envs/fetch.yaml'
    threads: get_threads('fetch_taxdump', 1)
    log:
        'logs/fetch_taxdump.log'
    benchmark:
        'logs/fetch_taxdump.benchmark.txt'
    resources:
        download=1,
        threads=get_threads('fetch_taxdump', 1)
    shell:
        '(wget -q -O taxdump.tar.gz ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz \
        && tar xzvf taxdump.tar.gz \
        && mv *.dmp {params.outdir} \
        && rm taxdump.tar.gz) 2> {log}'
