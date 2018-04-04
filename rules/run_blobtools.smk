rule blobtools_create:
    """
    Use BlobTools create to generate a BlobDB using the ordered similarity
    searches and coverage files.
    """
    input:
        assembly='{assembly}.fasta',
        dbs=list_similarity_results(config),
        coverage=expand('{{assembly}}.{sra}.bam.cov',sra=list_sra_accessions(reads))
    output:
        '{assembly}.blobDB.json'
    params:
        dbs=expand('-t {db}',db=list_similarity_results(config)),
        coverage=lambda wc: expand('-c '+wc.assembly+'.{sra}.bam.cov',sra=list_sra_accessions(reads)),
        taxrule=config['similarity']['taxrule'] if 'taxrule' in config['similarity'] else 'bestsumorder',
        assembly=lambda wc: wc.assembly
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
            -o {params.assembly}'
