rule run_windowmasker:
    """
    Run windowmasker to mask repeats in the assembly.
    """
    input:
        fasta = temp("{assembly}.fasta")
    output:
        counts = "{assembly}.windowmasker.counts",
        masked = "{assembly}.windowmasker.fasta"
    threads: 1
    log:
        "logs/{assembly}/run_windowmasker.log"
    benchmark:
        "logs/{assembly}/run_windowmasker.benchmark.txt"
    shell:
        """windowmasker -in {input.fasta} \
                      -infmt fasta \
                      -mk_counts \
                      -sformat obinary \
                      -out {output.counts} 2> {log} \
        && windowmasker -in {input} \
                        -infmt fasta \
                        -ustat {output.counts} \
                        -dust T \
                        -outfmt fasta \
                        -out {output.masked} 2>> {log}"""
