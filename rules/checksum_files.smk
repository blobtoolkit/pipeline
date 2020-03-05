rule checksum_files:
    """
    Calculate SHA1 checksum for all files in dataset.
    """
    input:
        '{assembly}/summary.json'
    output:
        '{assembly}/CHECKSUM'
    params:
        assembly = lambda wc: wc.assembly
    threads: get_threads('checksum_files', 1)
    log:
        'logs/{assembly}/checksum_files.log'
    benchmark:
        'logs/{assembly}/checksum_files.benchmark.txt'
    resources:
        threads = get_threads('checksum_files', 1)
    shell:
        '(find {params.assembly}/ -type f -exec sha1sum {{}} \';\' \
        | sort -k 2 \
        | sed \'s:{params.assembly}/::\' > {params.assembly}/CHECKSUM) 2> {log}'
