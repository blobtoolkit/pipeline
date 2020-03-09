rule subsample_fastq:
    """
    Subsample large fastq files to reduce mapping time.
    """
    input:
        lambda wc: "%s%s.gz" % (wc.sra, wc.suff.replace('.subsampled', ''))
    output:
        temp("{sra}{suff}.gz")
    params:
        cmd = lambda wc: generate_subsample_command(wc.sra, reads)
    wildcard_constraints:
        sra = r'[a-zA-Z0-9]+',
        suff = r'[\._\d\w]*\.subsampled.fast[aq]'
    conda:
        '../envs/pyblast.yaml'
    threads: get_threads('subsample_fastq', 1)
    log:
        "logs/%s/subsample_fastq/{sra}{suff}.log" % config['assembly']['prefix']
    benchmark:
        "logs/%s/subsample_fastq/{sra}{suff}.benchmark.txt" % config['assembly']['prefix']
    resources:
        download = 1,
        threads = get_threads('subsample_fastq', 1)
    shell:
        '({params.cmd[0]} {input} {params.cmd[1]} {output}) 2> {log}'
