rule map_reads:
    """
    Run bwa/minimap2 with settings appropriate to sequencing platform
    """
    input:
        fastq = lambda wc: list_read_files(wc.sra, reads, True),
        fasta = '{assembly}.fasta',
        index = expand('{{assembly}}.fasta.{suffix}', suffix=BWA_INDEX)
    output:
        '{assembly}.{sra}.bam' if keep else temp('{assembly}.{sra}.bam')
    params:
        cmd = lambda wc: generate_mapping_command(wc.sra, reads)
    conda:
        '../envs/bwa.yaml'
    threads: get_threads('map_reads', maxcore)
    log:
        'logs/{assembly}/map_reads/{sra}.log'
    benchmark:
        'logs/{assembly}/map_reads/{sra}.benchmark.txt'
    resources:
        threads = get_threads('map_reads', maxcore)
    shell:
        '({params.cmd} -t {threads} {input.fasta} {input.fastq} | \
        samtools view -h -T {input.fasta} - | \
        samtools sort -@{threads} -O BAM -o {output} -) 2> {log}'
