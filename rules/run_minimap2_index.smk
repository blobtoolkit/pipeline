rule run_minimap2_index:
    """
    Run minimap2 reference genome indexing.
    """
    input:
        config["assembly"]["file"]
    output:
        "{assembly}.{tuning}.mmi"
    params:
        tuning = lambda wc: wc.tuning,
        assembly = lambda wc: wc.assembly,
        span = config["assembly"]["span"]
    threads: 3
    log:
        "logs/{assembly}/run_minimap2_index/{tuning}.log"
    benchmark:
        "logs/{assembly}/run_minimap2_index/{tuning}.benchmark.txt"
    shell:
        """(minimap2 -t {threads} -x {params.tuning} -I {params.span} -d {output} {input}) &> {log}"""