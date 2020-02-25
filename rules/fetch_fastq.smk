rule fetch_fastq:
    """
    Fetch fastq file from EBI using aria2.
    """
    output:
        '{sra}{suff}.gz' if keep else temp('{sra}{suff}.gz')
    params:
        url = lambda wc: prepare_ebi_sra_url(wc.sra, "%s%s.gz" % (wc.sra, wc.suff))
    wildcard_constraints:
        sra = r'\wRR\d+',
        suff = r'[_\d\w]*\.fastq'
    conda:
        '../envs/fetch.yaml'
    threads: get_threads('fetch_fastq', 1)
    log:
        "logs/%s/fetch_fastq/{sra}{suff}.log" % config['assembly']['prefix']
    benchmark:
        "logs/%s/fetch_fastq/{sra}{suff}.benchmark.txt" % config['assembly']['prefix']
    resources:
        download = 1,
        threads = get_threads('fetch_fastq', 1)
    shell:
        'aria2c -c \
        --max-connection-per-server=8 \
        --min-split-size=1M \
        -s 8 \
        -l {log} \
        --log-level=notice \
        --show-console-readout=false \
        --console-log-level=error \
        {params.url}'
