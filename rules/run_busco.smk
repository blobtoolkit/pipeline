rule run_busco:
    """
    Run BUSCO on a lineage.
    """
    input:
        fasta='{assembly}.fasta',
        busco=config['busco']['lineage_dir']+'/{lineage}.tar.gz',
    output:
        full='{assembly}_{lineage}.tsv', # if keep else temp('{assembly}_{lineage}.tsv'),
        short='{assembly}_{lineage}_short_summary.txt' # if keep else temp('short_summary_{assembly}_{lineage}.txt')
    params:
        lineage=lambda wc: wc.lineage,
        lineage_dir=config['busco']['lineage_dir'],
        assembly=lambda wc: wc.assembly,
        outdir=lambda wc: "run_%s_%s" % (wc.assembly, wc.lineage)
    wildcard_constraints:
        lineage='\w+_odb9'
    conda:
        '../envs/busco.yaml'
    threads: 16
    resources:
        threads=16
    shell:
        'run_busco \
            -f \
            -i {input.fasta} \
            -o {params.assembly}_{params.lineage} \
            -l {params.lineage_dir}/{params.lineage} \
            -m geno \
            -c {threads} && \
        mv {params.outdir}/full_table_{params.assembly}_{params.lineage}.tsv {output.full} && \
        mv {params.outdir}/short_summary_{params.assembly}_{params.lineage}.txt {output.short} && \
        rm -r {params.outdir} && exit 0 || \
        rm -r {params.outdir} && rm -r {params.assembly}_{params.lineage} exit 1'
