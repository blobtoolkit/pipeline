BWA_INDEX = ['amb','ann','bwt','pac','sa']

rule bwa_index:
    """
    Index an assembly FASTA file for use with BWA
    """
    input:
        '{assembly}.fna'
    output:
        temp(expand('{{assembly}}.fna.{suffix}',suffix=BWA_INDEX))
    conda:
         '../envs/bwa.yaml'
    threads: 1
    resources:
        threads=1
    shell:
        'bwa index -a bwtsw {input}'

rule interleaved_bwa_mem:
    """
    Run bwa mem with an interleaved FASTQ file
    """
    input:
        fastq='{sample}.interleaved.fastq',
        fna='{assembly}.fna',
        index=expand('{{assembly}}.fna.{suffix}',suffix=BWA_INDEX)
    output:
        temp('{assembly}.{sample}.bam')
    conda:
         '../envs/bwa.yaml'
    threads: 32
    resources:
        threads=32
    shell:
        'bwa mem -M -t {threads} -p {input.fna} {input.fastq} | \
        samtools sort -O BAM -o {output} -'

rule blobtools_map2cov:
    """
    Use BlobTools to convert a mapping file to a coverage file
    """
    input:
        bam='{assembly}.{sample}.bam',
        fna='{assembly}.fna'
    output:
        temp('{assembly}.{sample}.bam.cov')
    conda:
         '../envs/blobtools.yaml'
    threads: 1
    resources:
        threads=1
    shell:
        'blobtools map2cov \
            -i {input.fna} \
            -b {input.bam}'
