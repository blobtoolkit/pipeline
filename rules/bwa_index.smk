rule bwa_index:
    """
    Index an assembly FASTA file for use with BWA
    """
    input:
        '{assembly}.fasta'
    output:
        temp(expand('{{assembly}}.fasta.{suffix}', suffix=BWA_INDEX))
    conda:
        '../envs/bwa.yaml'
    threads: get_threads('bwa_index', 1)
    log:
        'logs/{assembly}/bwa_index.log'
    benchmark:
        'logs/{assembly}/bwa_index.benchmark.txt'
    resources:
        threads = get_threads('bwa_index', 1)
    shell:
        'bwa index -a bwtsw {input} > {log} 2>&1'
