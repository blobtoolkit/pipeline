rule blobtoolkit_create:
    """
    Use BlobTools2 create to generate a BlobDir using an assembly fasta file and metadata.
    """
    input:
        fasta = '{assembly}.fasta',
        yaml = '{assembly}.meta.yaml',
        lineages = "%s/taxidlineage.dmp" % config['settings']['taxonomy']
    output:
        "{assembly}%s/identifiers.json" % rev
    params:
        taxdump = taxdump_dir,
        id = lambda wc: "%s%s" % (wc.assembly, rev),
        taxid = config['taxon']['taxid'],
        path = config['settings']['blobtools2_path']
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('blobtoolkit_create', 1)
    log:
        'logs/{assembly}/blobtoolkit_create.log'
    benchmark:
        'logs/{assembly}/blobtoolkit_create.benchmark.txt'
    resources:
        btk = 1
    shell:
        'PATH={params.path}:$PATH && \
        blobtools replace \
            --fasta {input.fasta} \
            --meta {input.yaml} \
            --taxdump "{params.taxdump}" \
            --taxid {params.taxid} \
            {params.id} > {log} 2>&1'
