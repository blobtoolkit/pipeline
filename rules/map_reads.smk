rule bwa_index:
    """
    Index an assembly FASTA file for use with BWA
    """
    input:
        '{assembly}.fasta'
    output:
        temp(expand('{{assembly}}.fasta.{suffix}',suffix=BWA_INDEX))
    conda:
        '../envs/bwa.yaml'
    threads: get_threads('bwa_index', 1)
    log:
        lambda wc: "logs/%s/bwa_index.log" % (wc.assembly)
    benchmark:
        'logs/{assembly}/bwa_index.benchmark.txt'
    resources:
        threads=get_threads('bwa_index', 1)
    shell:
        'bwa index -a bwtsw {input} > {log} 2>&1'

rule map_reads:
    """
    Run bwa/minimap2 with settings appropriate to sequencing platform
    """
    input:
        fastq=lambda wc: list_read_files(wc.sra,reads,True),
        fasta='{assembly}.fasta',
        index=expand('{{assembly}}.fasta.{suffix}',suffix=BWA_INDEX)
    output:
        '{assembly}.{sra}.bam' if keep else temp('{assembly}.{sra}.bam')
    params:
        cmd = lambda wc: generate_mapping_command(wc.sra,reads)
    conda:
        '../envs/bwa.yaml'
    threads: get_threads('map_reads', maxcore)
    log:
        lambda wc: "logs/%s/map_reads/%s.log" % (wc.assembly, wc.sra)
    benchmark:
        'logs/{assembly}/map_reads/{sra}.benchmark.txt'
    resources:
        threads=get_threads('map_reads', maxcore)
    shell:
        '({params.cmd} -t {threads} {input.fasta} {input.fastq} | \
        samtools view -h -T {input.fasta} - | \
        samtools sort -@{threads} -O BAM -o {output} -) 2> {log}'

rule bamtools_stats:
    """
    Run bamtools stats to generate summary statistics for each BAM file
    """
    input:
        bam='{assembly}.{sra}.bam'
    output:
        '{assembly}.{sra}.bam.stats'
    conda:
        '../envs/bwa.yaml'
    threads: get_threads('bamtools_stats', 1)
    log:
        lambda wc: "logs/%s/bamtools_stats/%s.log" % (wc.assembly, wc.sra)
    benchmark:
        'logs/{assembly}/bamtools_stats/{sra}.benchmark.txt'
    resources:
        threads=get_threads('bamtools_stats', 1)
    shell:
        'bamtools stats \
            -in {input.bam} \
            -insert \
            > {output} 2> {log}'
