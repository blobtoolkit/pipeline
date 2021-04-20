rule run_mosdepth:
    """
    Run Mosdepth to get coverage depth.
    """
    input:
        bam = "%s/{assembly}.{sra}.bam" % minimap_path,
        csi = "%s/{assembly}.{sra}.bam.csi" % minimap_path,
        bed = "{assembly}.stats.mask.bed"
    output:
        bed = "{assembly}.stats.{sra}_cov.bed",
        bed_windows = "{assembly}.stats.{sra}_cov_windows.bed"
    params:
        prefix = lambda wc: "%s.%s" % (wc.assembly, wc.sra)
    threads: 4
    log:
        "logs/{assembly}/run_mosdepth/{sra}.log"
    benchmark:
        "logs/{assembly}/run_mosdepth/{sra}.benchmark.txt"
    shell:
        """(mosdepth -n -x -b {input.bed} \
                 -t 4 {params.prefix} {input.bam} && \
        > {output.bed} && \
        zcat {params.prefix}.regions.bed.gz | grep full | perl -lne '@x=split/\s+/;$x[3]=".";print join("\t",@x)' >> {output.bed} && \
        > {output.bed_windows} && \
        zcat {params.prefix}.regions.bed.gz | grep window | perl -lne '@x=split/\s+/;$x[3]=".";print join("\t",@x)' >> {output.bed_windows} \
        ) 2> {log}"""