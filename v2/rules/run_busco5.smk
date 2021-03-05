rule run_busco_v5:
    """
    Run BUSCO on a lineage.
    """
    input:
        fasta = "{assembly}.fasta",
    output:
        full = "{assembly}.busco.{lineage}.tsv",
        short = "{assembly}.busco.{lineage}.txt"
    params:
        lineage_dir = lambda wc: "%s/%s" % (config["busco"]["lineage_dir"], wc.lineage),
        lineage = lambda wc: wc.lineage,
        assembly = lambda wc: wc.assembly,
        outdir = lambda wc: "%s_%s" % (wc.assembly, wc.lineage)
    threads: 30
    log:
        "logs/{assembly}/run_busco/{lineage}.log"
    benchmark:
        "logs/{assembly}/run_busco/{lineage}.benchmark.txt"
    shell:
        """busco \
            -f \
            -i {input.fasta} \
            -o {params.assembly}_{params.lineage} \
            -l {params.lineage_dir} \
            -m geno \
            -c {threads} > {log} 2>&1 && \
        mv {params.outdir}/run_{params.lineage}/full_table.tsv {output.full} && \
        mv {params.outdir}/run_{params.lineage}/short_summary.txt {output.short} && \
        rm -rf {params.outdir} && exit 0 || \
        rm -rf {params.outdir} && exit 1"""