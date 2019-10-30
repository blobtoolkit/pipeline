rule fetch_assembly:
    """
    Fetch a remote assembly from EBI or NCBI.
    """
    output:
        fa='{assembly}.fasta'
    params:
        url=lambda wc: prepare_ncbi_assembly_url(config['assembly']['accession'],config['assembly']['alias'])
    wildcard_constraints:
        assembly='\w+'
    conda:
         '../envs/fetch.yaml'
    threads: 1
    log:
        lambda wc: "logs/%s/fetch_assembly.log" % (wc.assembly)
    benchmark:
        'logs/{assembly}/fetch_assembly.benchmark.txt'
    resources:
        download=1,
        threads=1
    shell:
        '(curl -s {params.url} | \
        pigz -d > {output.fa}) 2> {log}'

rule fetch_fastq:
    """
    Fetch fastq file from EBI using aria2.
    """
    output:
        '{sra}{suff}.gz' if keep else temp('{sra}{suff}.gz')
    params:
        url = lambda wc: prepare_ebi_sra_url(wc.sra,"%s%s.gz" % (wc.sra, wc.suff))
    wildcard_constraints:
        sra='\wRR\d+',
        suff='[_\d\w]*\.fastq'
    conda:
         '../envs/fetch.yaml'
    threads: 1
    log:
        lambda wc: "logs/%s/fetch_fastq/%s%s.log" % (config['assembly']['prefix'], wc.sra, wc.suff)
    benchmark:
        "logs/%s/fetch_fastq/{sra}{suff}.benchmark.txt" % config['assembly']['prefix']
    resources:
        download=1,
        threads=1
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

rule subsample_fastq:
    """
    Subsample large fastq files to reduce mapping time.
    """
    input:
        lambda wc: "%s%s.gz" % (wc.sra,wc.suff.replace('.subsampled',''))
    output:
        temp("{sra}{suff}.gz")
    params:
        cmd = lambda wc: generate_subsample_command(wc.sra,reads)
    wildcard_constraints:
        sra='\wRR\d+',
        suff='[_\dsubread]*\.subsampled.fastq'
    conda:
         '../envs/pyblast.yaml'
    threads: 1
    log:
        lambda wc: "logs/%s/subsample_fastq/%s%s.log" % (config['assembly']['prefix'], wc.sra, wc.suff)
    benchmark:
        "logs/%s/subsample_fastq/{sra}{suff}.benchmark.txt" % config['assembly']['prefix']
    resources:
        download=1,
        threads=1
    shell:
        '({params.cmd[0]} {input} {params.cmd[1]} {output}) 2> {log}'
