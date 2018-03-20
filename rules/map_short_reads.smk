rule bwa_index:
    """
    Index an assembly FASTA file for use with BWA
    """
    input:
        '{assembly}.fna'
    output:
        '{assembly}.fna.bwt'
    conda:
         '../envs/bwa.yaml'
    threads: 1
    shell:
        '{ENV} \
        bwa index -a bwtsw {input}'

rule interleaved_bwa_mem:
    """
    Run bwa mem with an interleaved FASTQ file
    """
    input:
        fastq='{sample}.interleaved.fastq',
        fna='{assembly}.fna',
        index='{assembly}.fna.bwt'
    output:
        '{assembly}.{sample}.bam'
    conda:
         '../envs/bwa.yaml'
    threads: 32
    shell:
        '{ENV} \
        bwa mem -M -t {threads} -p {input.fna} {input.fastq} | \
        samtools sort -O BAM -o {output} -'

rule blobtools_map2cov:
    """
    Use BlobTools to convert a mapping file to a coverage file
    """
    input:
        bam='{assembly}.{sample}.bam',
        fna='{assembly}.fna'
    output:
        '{assembly}.{sample}.bam.cov'
    conda:
         '../envs/blobtools.yaml'
    threads: 1
    shell:
        '{ENV} \
        PATH=' + config['settings']['blobtools'] + ':$PATH && \
        blobtools map2cov \
            -i {input.fna} \
            -b {input.bam}'
