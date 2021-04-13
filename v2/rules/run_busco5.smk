rule run_busco_v5:
    """
    Run BUSCO on a lineage.
    """
    input:
        fasta = "{assembly}.fasta",
    output:
        full = "{assembly}.busco.{lineage}/full_table.tsv.gz",
    params:
        lineage_dir = lambda wc: "%s/%s" % (config["busco"]["lineage_dir"], wc.lineage),
        lineage = lambda wc: wc.lineage,
        assembly = lambda wc: wc.assembly,
        outdir = lambda wc: "%s_%s" % (wc.assembly, wc.lineage),
        buscodir = lambda wc: "%s.busco.%s" % (wc.assembly, wc.lineage),
        missing = lambda wc: "%s.busco.%s/missing_busco_list.tsv.gz" % (wc.assembly, wc.lineage),
        short = lambda wc: "%s.busco.%s/short_summary.txt.gz" % (wc.assembly, wc.lineage),
        sequences = lambda wc: "%s.busco.%s/busco_sequences.tar.gz" % (wc.assembly, wc.lineage)
    threads: 30
    log:
        "logs/{assembly}/run_busco/{lineage}.log"
    benchmark:
        "logs/{assembly}/run_busco/{lineage}.benchmark.txt"
    shell:
        """mkdir -p {params.buscodir} && \
        busco \
            -f \
            -i {input.fasta} \
            -o {params.assembly}_{params.lineage} \
            -l {params.lineage_dir} \
            -m geno \
            -c {threads} > {log} 2>&1 && \
        gzip -c {params.outdir}/run_{params.lineage}/full_table.tsv > {output.full} && \
        gzip -c {params.outdir}/run_{params.lineage}/short_summary.txt > {params.short} && \
        gzip -c {params.outdir}/run_{params.lineage}/missing_busco_list.tsv > {params.missing} && \
        tar czf {params.sequences} -C {params.outdir}/run_{params.lineage} busco_sequences && \
        rm -rf {params.outdir} && exit 0 || \
        if grep -q "Prodigal did not recognize any genes" {log}; then \
            touch {output.full} && \
            rm -rf {params.outdir} && exit 0; \
        fi; \
        rm -rf {params.outdir} && exit 1"""

# ERROR:  Prodigal did not recognize any genes matching the dataset archaea_odb10 in the input file. If this is unexpected, check your input file and your installation of Prodigal
