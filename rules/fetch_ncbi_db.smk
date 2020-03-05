rule fetch_ncbi_db:
    """
    Fetch formatted NCBI BLAST databases
    """
    output:
        "%s/{name}.{suffix}" % similarity['blastdb']['local']
    wildcard_constraints:
        # NB: the path to the local copy of the file must contain the string 'ncbi'
        path = '.+ncbi.+',
        suffix = '[np]al'
    params:
        ftp_url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/',
        ftp_dir = lambda wc: "%s" % 'v5/' if wc.name.endswith('_v5') else '',
        name = lambda wc: wc.name,
        path = ncbi_dir
    conda:
        '../envs/fetch.yaml'
    threads: get_threads('fetch_ncbi_db', 1)
    # log:
    #     lambda wc: "logs/fetch_ncbi_db/%s.%s.log" % (wc.name, wc.suffix)
    log:
        'logs/fetch_ncbi_db/{name}.{suffix}.log'
    benchmark:
        'logs/fetch_ncbi_db/{name}.{suffix}.benchmark.txt'
    resources:
        download = 1,
        threads = get_threads('fetch_ncbi_db', 1)
    shell:
        '(wget "{params.ftp_url}{params.ftp_dir}{params.name}.??.tar.gz" -P {params.path}/ && \
        for file in {params.path}/*.tar.gz; \
            do tar xf $file -C {params.path} && rm $file; \
        done) 2> {log}'
