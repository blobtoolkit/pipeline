rule fetch_fastq:
    """
    Fetch fastq file from EBI using aria2.
    """
    output:
        '{sra}{suff}.gz' if keep else temp('{sra}{suff}.gz')
    params:
        url = lambda wc: prepare_ebi_sra_url(wc.sra, "%s%s.gz" % (wc.sra, wc.suff))
    wildcard_constraints:
        sra = r'[a-zA-Z0-9]+',
        suff = r'(?!.*\b(?:subsampled)\b)[\._\d\w]*\.fast[aq]'
    conda:
        '../envs/fetch.yaml'
    threads: get_threads('fetch_fastq', 1)
    log:
        "logs/%s/fetch_fastq/{sra}{suff}.log" % config['assembly']['prefix']
    benchmark:
        "logs/%s/fetch_fastq/{sra}{suff}.benchmark.txt" % config['assembly']['prefix']
    resources:
        download = 1
    shell:
        'if [[ "{params.url}" == *tp://* ]]; then \
            aria2c -c \
            --max-connection-per-server=8 \
            --min-split-size=1M \
            -s 8 \
            -l {log} \
            --log-level=notice \
            --show-console-readout=false \
            --console-log-level=error \
            {params.url}; \
        else \
            cp -n {params.url} ./ 2> {log}; \
        fi'
