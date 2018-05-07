rule blobtools_create:
    """
    Use BlobTools create to generate a BlobDB using the ordered similarity
    searches and coverage files.
    """
    input:
        assembly='{assembly}.fasta',
        dbs=list_similarity_results(config),
        coverage=expand('{{assembly}}.{sra}.bam.cov',sra=list_sra_accessions(reads)),
        covsum=lambda wc:platform_cov_files(reads,wc.assembly)
    output:
        '{assembly}.blobDB.json'
    params:
        dbs=expand('-t {db}',db=list_similarity_results(config)),
        coverage=lambda wc: expand('-c '+wc.assembly+'.{sra}.bam.cov',sra=list_sra_accessions(reads)),
        taxrule=config['similarity']['taxrule'] if 'taxrule' in config['similarity'] else 'bestsumorder',
        assembly=lambda wc: wc.assembly,
        covsum=lambda wc:list(map(lambda file: '-c '+file, platform_cov_files(reads,wc.assembly)))
    conda:
        '../envs/blobtools.yaml'
    threads: 1
    resources:
        threads=1
    shell:
        'blobtools create \
            -i {input.assembly} \
            {params.dbs} \
            -x "{params.taxrule}" \
            {params.coverage} \
            {params.covsum} \
            -o {params.assembly}'
