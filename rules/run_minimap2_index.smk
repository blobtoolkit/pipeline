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
        assembly = lambda wc: wc.assembly
    threads: 3
    log:
        "logs/{assembly}/run_minimap2_index/{tuning}.log"
    benchmark:
        "logs/{assembly}/run_minimap2_index/{tuning}.benchmark.txt"
    shell:
        """(minimap2 -t {threads} -x {params.tuning} -d {output} {input}) &> {log}"""