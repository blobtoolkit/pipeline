import os

rule generate_metadata:
    input:
        stats=expand('{{assembly}}.{sra}.bam.stats',sra=list_sra_accessions(reads)),
        yaml='{assembly}.yaml'
    output:
        '{assembly}.meta.yaml'
    conda:
        '../envs/py3.yaml'
    params:
        gitdir=os.path.dirname(os.path.abspath(workflow.snakefile))+'/.git'
    threads: 1
    resources:
        threads=1
    script:
        '../scripts/generate_metadata.py'


rule blobtoolkit_create:
    """
    Use BlobTools2 create to generate a BlobDir using the ordered similarity
    searches and assembly files.
    """
    input:
        fasta='{assembly}.fasta',
        yaml='{assembly}.meta.yaml',
        dbs=list_similarity_results(config),
        lineages="%s/taxidlineage.dmp" % config['settings']['taxonomy']
    output:
        '{assembly}/meta.json'
    params:
        taxrule=config['similarity']['taxrule'] if 'taxrule' in config['similarity'] else 'bestsumorder',
        taxdump=config['settings']['taxonomy'],
        assembly=lambda wc: wc.assembly,
        path=config['settings']['blobtools2_path'],
        taxid=config['taxon']['taxid'],
        dbs='.raw --hits '.join(list_similarity_results(config))
    conda:
        '../envs/blobtools2.yaml'
    threads: 1
    resources:
        threads=1
    shell:
        '{params.path}/blobtools create \
            --fasta {input.fasta} \
            --meta {input.yaml} \
            --hits {params.dbs} \
            --taxrule "{params.taxrule}" \
            --taxdump "{params.taxdump}" \
            --taxid {params.taxid} \
            {params.assembly}'

rule blobtoolkit_add_cov:
    """
    Use BlobTools2 add to add coverage to a BlobDir from BAM files.
    """
    input:
        meta="%s/meta.json" % config['assembly']['prefix'],
        bam=expand("%s.{sra}.bam" % asm, sra=list_sra_accessions(reads))
    output:
        expand("%s/{sra}_cov.json" % config['assembly']['prefix'],sra=list_sra_accessions(reads))
    params:
        assembly=config['assembly']['prefix'],
        path=config['settings']['blobtools2_path'],
        covs=lambda wc: ' --cov '.join(["%s.%s.bam" % (config['assembly']['prefix'], sra) for sra in list_sra_accessions(reads)])
    conda:
        '../envs/blobtools2.yaml'
    threads: 1
    resources:
        threads=1
    shell:
        '{params.path}/blobtools add \
            --cov {params.covs} \
            {params.assembly}'


