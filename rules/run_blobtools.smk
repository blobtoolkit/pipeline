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
    log:
      lambda wc: "logs/%s/generate_metadata.log" % (wc.assembly)
    resources:
        threads=1
    script:
        '../scripts/generate_metadata.py'


rule blobtoolkit_create:
    """
    Use BlobTools2 create to generate a BlobDir using an assembly fasta file and metadata.
    """
    input:
        fasta='{assembly}.fasta',
        yaml='{assembly}.meta.yaml',
        lineages="%s/taxidlineage.dmp" % config['settings']['taxonomy']
    output:
        '{assembly}/identifiers.json'
    params:
        taxdump=config['settings']['taxonomy'],
        assembly=lambda wc: wc.assembly,
        path=config['settings']['blobtools2_path'],
        taxid=config['taxon']['taxid'],
    conda:
        '../envs/blobtools2.yaml'
    threads: 1
    log:
      lambda wc: "logs/%s/blobtoolkit_create.log" % (wc.assembly)
    resources:
        threads=1,
        btk=1
    shell:
        '{params.path}/blobtools create \
            --fasta {input.fasta} \
            --meta {input.yaml} \
            --taxdump "{params.taxdump}" \
            --taxid {params.taxid} \
            {params.assembly} > {log} 2>&1'

rule blobtoolkit_add_hits:
    """
    Add ordered similarity search results to a BlobDir.
    """
    input:
        meta="%s/identifiers.json" % config['assembly']['prefix'],
        dbs=list_similarity_results(config),
        lineages="%s/taxidlineage.dmp" % config['settings']['taxonomy']
    output:
        "{assembly}/%s_phylum_positions.json" % config['similarity']['taxrule']
    params:
        taxrule=config['similarity']['taxrule'] if 'taxrule' in config['similarity'] else 'bestsumorder',
        taxdump=config['settings']['taxonomy'],
        assembly=lambda wc: wc.assembly,
        path=config['settings']['blobtools2_path'],
        dbs='.raw --hits '.join(list_similarity_results(config))
    conda:
        '../envs/blobtools2.yaml'
    threads: 1
    log:
      lambda wc: "logs/%s/blobtoolkit_add_hits.log" % (wc.assembly)
    resources:
        threads=1,
        btk=1
    shell:
        '{params.path}/blobtools add \
            --hits {params.dbs} \
            --taxrule "{params.taxrule}" \
            --taxdump "{params.taxdump}" \
            {params.assembly} > {log} 2>&1'

rule blobtoolkit_add_cov:
    """
    Use BlobTools2 add to add coverage to a BlobDir from BAM files.
    """
    input:
        meta="%s/identifiers.json" % config['assembly']['prefix'],
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
    log:
      lambda wc: "logs/%s/blobtoolkit_add_cov.log" % (config['assembly']['prefix'])
    resources:
        threads=1,
        btk=1
    shell:
        '{params.path}/blobtools add \
            --cov {params.covs} \
            --threads {threads} \
            {params.assembly} > {log} 2>&1'


rule blobtoolkit_add_busco:
    """
    import BUSCO results into BlobDir.
    """
    input:
        meta="%s/identifiers.json" % config['assembly']['prefix'],
        tsv=expand("%s_{lineage}.tsv" % config['assembly']['prefix'],lineage=config['busco']['lineages'])
    output:
        expand("%s/{lineage}_busco.json" % config['assembly']['prefix'],lineage=config['busco']['lineages'])
    params:
        assembly=config['assembly']['prefix'],
        path=config['settings']['blobtools2_path'],
        busco=' --busco '.join(["%s_%s.tsv" % (config['assembly']['prefix'],lineage) for lineage in config['busco']['lineages']])
    conda:
        '../envs/blobtools2.yaml'
    threads: 1
    log:
      lambda wc: "logs/%s/blobtoolkit_add_busco.log" % (config['assembly']['prefix'])
    resources:
        threads=1,
        btk=1
    shell:
        '{params.path}/blobtools replace \
            --busco {params.busco} \
            {params.assembly} > {log} 2>&1'
