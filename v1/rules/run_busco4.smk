rule run_busco_v4:
    """
    Run BUSCO on a lineage.
    """
    input:
        fasta = '{assembly}.fasta',
    output:
        full = '{assembly}.busco.{lineage}.tsv',
        short = '{assembly}.busco.{lineage}.txt'
    params:
        lineage = lambda wc: wc.lineage,
        assembly = lambda wc: wc.assembly,
        outdir = lambda wc: "%s_%s" % (wc.assembly, wc.lineage)
    wildcard_constraints:
        lineage = r'\w+_odb10'
    conda:
        '../envs/busco4.yaml'
    threads: get_threads('run_busco', multicore)
    log:
        'logs/{assembly}/run_busco/{lineage}.log'
    benchmark:
        'logs/{assembly}/run_busco/{lineage}.benchmark.txt'
    resources:
        threads = get_threads('run_busco', multicore)
    shell:
        'busco \
            -f \
            -i {input.fasta} \
            -o {params.assembly}_{params.lineage} \
            -l {params.lineage} \
            -m geno \
            -c {threads} > {log} 2>&1 && \
        mv {params.outdir}/run_{params.lineage}/full_table.tsv {output.full} && \
        mv {params.outdir}/run_{params.lineage}/short_summary.txt {output.short} && \
        rm -rf {params.outdir} && exit 0 || \
        rm -rf {params.outdir} && exit 1'
