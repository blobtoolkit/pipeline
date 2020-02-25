rule extract_nohit_sequences:
    """
    Run seqtk to extract nohit sequences from assembly.
    """
    input:
        fasta = '{assembly}.fasta',
        nohit = '{assembly}.blastn.{name}.root.{root}{masked}.nohit'
    output:
        '{assembly}.blastn.{name}.root.{root}{masked}.fasta.nohit' if keep else temp('{assembly}.blastn.{name}.root.{root}{masked}.fasta.nohit')
    conda:
        '../envs/pyblast.yaml'
    threads: get_threads('extract_nohit_sequences', 1)
    log:
        'logs/{assembly}/extract_nohit_sequences/{name}.root.{root}{masked}.log'
    benchmark:
        'logs/{assembly}/extract_nohit_sequences/{name}.root.{root}{masked}.benchmark.txt'
    resources:
        threads = get_threads('extract_nohit_sequences', 1)
    shell:
        'seqtk subseq {input.fasta} {input.nohit} > {output} 2> {log}'
