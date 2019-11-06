rule run_busco:
    """
    Run BUSCO on a lineage.
    """
    input:
        fasta='{assembly}.fasta',
        busco=config['busco']['lineage_dir']+'/{lineage}.tar.gz',
    output:
        full='{assembly}.{lineage}.tsv',
        short='{assembly}.{lineage}.txt'
    params:
        lineage=lambda wc: wc.lineage,
        lineage_dir=config['busco']['lineage_dir'],
        assembly=lambda wc: wc.assembly,
        outdir=lambda wc: "run_%s_%s" % (wc.assembly, wc.lineage)
    wildcard_constraints:
        lineage='\w+_odb9'
    conda:
        '../envs/busco.yaml'
    threads: get_threads('run_busco', multicore)
    log:
        lambda wc: "logs/%s/run_busco/%s.log" % (wc.assembly, wc.lineage)
    benchmark:
        'logs/{assembly}/run_busco/{lineage}.benchmark.txt'
    resources:
        threads=get_threads('run_busco', multicore)
    shell:
        'run_busco \
            -f \
            -i {input.fasta} \
            -o {params.assembly}_{params.lineage} \
            -l {params.lineage_dir}/{params.lineage} \
            -m geno \
            -c {threads} > {log} 2>&1 && \
        mv {params.outdir}/full_table_{params.assembly}_{params.lineage}.tsv {output.full} && \
        mv {params.outdir}/short_summary_{params.assembly}_{params.lineage}.txt {output.short} && \
        rm -rf {params.outdir} && rm -rf {params.assembly}_{params.lineage} && exit 0 || \
        rm -rf {params.outdir} && rm -rf {params.assembly}_{params.lineage} exit 1'
