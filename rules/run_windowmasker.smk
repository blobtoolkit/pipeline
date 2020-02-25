rule run_windowmasker:
    """
    Run windowmasker to mask repeats in the assembly.
    """
    input:
        '{assembly}.fasta'
    output:
        counts = '{assembly}.windowmasker.counts' if keep else temp('{assembly}.windowmasker.counts'),
        masked = '{assembly}.windowmasker.fasta' if keep else temp('{assembly}.windowmasker.fasta')
    conda:
        '../envs/pyblast.yaml'
    threads: get_threads('run_windowmasker', 1)
    log:
        'logs/{assembly}/run_windowmasker.log'
    benchmark:
        'logs/{assembly}/run_windowmasker.benchmark.txt'
    resources:
        threads = get_threads('run_windowmasker', 1)
    shell:
        'windowmasker -in {input} \
                      -infmt fasta \
                      -mk_counts \
                      -sformat obinary \
                      -out {output.counts} 2> {log} \
        && windowmasker -in {input} \
                        -infmt fasta \
                        -ustat {output.counts} \
                        -dust T \
                        -outfmt fasta \
                        -out {output.masked} 2>> {log} '
