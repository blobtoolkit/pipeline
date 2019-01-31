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
        covsum=lambda wc:list(map(lambda file: '-c '+file, platform_cov_files(reads,wc.assembly))),
        path=config['settings']['blobtools_path']
    conda:
        '../envs/blobtools.yaml'
    threads: 1
    resources:
        threads=1
    shell:
        '{params.path}/blobtools create \
            -i {input.assembly} \
            {params.dbs} \
            -x "{params.taxrule}" \
            {params.coverage} \
            {params.covsum} \
            -o {params.assembly}'

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

rule blobtools_view:
    """
    Use BlobTools create to generate a BlobDB using the ordered similarity
    searches and coverage files.
    """
    input:
        blobDB='{assembly}.blobDB.json',
        yaml='{assembly}.meta.yaml'
    output:
        '{assembly}/meta.json'
    params:
        taxrule=config['similarity']['taxrule'] if 'taxrule' in config['similarity'] else 'bestsumorder',
        assembly=lambda wc: wc.assembly,
        path=config['settings']['blobtools_path']
    conda:
        '../envs/blobtools.yaml'
    threads: 1
    resources:
        threads=1
    shell:
        '{params.path}/blobtools view \
            -i {input.blobDB} \
            -x "{params.taxrule}" \
            --experimental {input.yaml}'

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
        dbs='--hits '.join(list_similarity_results(config))
    conda:
        '../envs/blobtools.yaml'
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
            --taxid {params.taxid}
            {params.assembly}'
