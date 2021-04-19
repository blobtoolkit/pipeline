rule generate_summary:
    """
    Use BlobTools2 filter to generate a dataset summary.
    """
    input:
        "{blobdir}/cumulative.png"
    output:
        "{blobdir}/summary.json"
    params:
        blobdir = lambda wc: wc.blobdir,
        taxrule = similarity_setting(config, "diamond_blastx", "taxrule")
    threads: 1
    log:
        "logs/{blobdir}/generate_summary.log"
    benchmark:
        "logs/{blobdir}/generate_summary.benchmark.txt"
    shell:
        """blobtools filter --summary {output} --taxrule {params.taxrule} {params.blobdir} 2> {log}"""
