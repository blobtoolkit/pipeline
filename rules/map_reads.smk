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
    threads: 1
    resources:
        threads=1
    shell:
        'bwa index -a bwtsw {input}'

rule map_reads:
    """
    Run bwa/minimap2 with settings appropriate to sequencing platform
    """
    input:
        fastq=lambda wc: list_read_files(wc.sra,reads,True),
        fasta='{assembly}.fasta',
        index=expand('{{assembly}}.fasta.{suffix}',suffix=BWA_INDEX)
    output:
        temp('{assembly}.{sra}.bam')
    params:
        cmd = lambda wc: generate_mapping_command(wc.sra,reads)
    conda:
         '../envs/bwa.yaml'
    threads: 32
    resources:
        threads=32
    shell:
        '{params.cmd} -t {threads} {input.fasta} {input.fastq} | \
        samtools sort -@{threads} -O BAM -o {output} -'

rule bamtools_stats:
    """
    Run bamtools stats to generate summary statistics for each BAM file
    """
    input:
        bam='{assembly}.{sra}.bam'
    output:
        '{assembly}.{sra}.bam.stats'
    conda:
         '../envs/blobtools.yaml'
    threads: 1
    resources:
        threads=1
    shell:
        'bamtools stats \
            -in {input.bam} \
            -insert \
            > {output}'

rule blobtools_map2cov:
    """
    Use BlobTools to convert a mapping file to a coverage file
    """
    input:
        bam='{assembly}.{sra}.bam',
        fasta='{assembly}.fasta'
    output:
        temp('{assembly}.{sra}.bam.cov')
    conda:
         '../envs/blobtools.yaml'
    threads: 1
    resources:
        threads=1
    shell:
        'blobtools map2cov \
            -i {input.fasta} \
            -b {input.bam}'

rule sum_coverage:
    """
    Sum coverage across a set of blobtools cov files
    """
    input:
        lambda wc: cov_files_by_platform(reads,wc.assembly,wc.platform)
    output:
        temp('{assembly}.{platform}.sum.cov')
    conda:
         '../envs/py3.yaml'
    threads: 1
    resources:
        threads=1
    script:
        '../scripts/sum_columns.py'
