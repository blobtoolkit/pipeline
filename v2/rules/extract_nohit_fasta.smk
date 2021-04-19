rule extract_nohit_fasta:
    """
    Extract sequences with no blastx hits into a separate file.
    """
    input:
        blastx = "%s/%s.diamond.reference_proteomes.out" % (diamond_path, config["assembly"]["prefix"]),
        fasta = "%s/{assembly}.windowmasker.fasta" % windowmasker_path
    output:
        "{assembly}.nohit.fasta"
    threads: 4
    log:
        "logs/{assembly}/extract_nohit_fasta.log"
    benchmark:
        "logs/{assembly}/extract_nohit_fasta.benchmark.txt"
    shell:
        """seqtk subseq {input.fasta} <(grep '<' {input.fasta} | \
            grep -v -f -F <(cut -f1 {input.blastx} | sort | uniq)) > {output} 2> {log}"""
