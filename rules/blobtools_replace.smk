import os

rule blobtoolkit_replace_hits:
    """
    Add ordered similarity search results to a BlobDir.
    """
    input:
        meta="%s%s/identifiers.json" % (config['assembly']['prefix'],vers),
        dbs=list_similarity_results(config),
        lineages="%s/taxidlineage.dmp" % (config['settings']['taxonomy'])
    output:
        "{assembly}%s/%s_phylum_positions.json" % (vers,config['similarity']['taxrule'])
    params:
        taxrule=config['similarity']['taxrule'] if 'taxrule' in config['similarity'] else 'bestsumorder',
        taxdump=config['settings']['taxonomy'],
        id=lambda wc: "%s%s" % (wc.assembly,vers),
        path=config['settings']['blobtools2_path'],
        dbs='.raw --hits '.join(list_similarity_results(config))
    conda:
        '../envs/blobtools2.yaml'
    threads: 1
    log:
      lambda wc: "logs/%s/blobtoolkit_replace_hits.log" % (wc.assembly)
    resources:
        threads=1,
        btk=1
    shell:
        '{params.path}/blobtools replace \
            --hits {params.dbs} \
            --taxrule "{params.taxrule}" \
            --taxdump "{params.taxdump}" \
            {params.id} > {log} 2>&1'

rule blobtoolkit_replace_cov:
    """
    Use BlobTools2 add to add coverage to a BlobDir from BAM files.
    """
    input:
        meta="%s%s/identifiers.json" % (config['assembly']['prefix'],vers),
        bam=expand("%s.{sra}.bam" % asm, sra=list_sra_accessions(reads))
    output:
        expand("%s%s/{sra}_cov.json" % (config['assembly']['prefix'],vers),sra=list_sra_accessions(reads))
    params:
        id="%s%s" % (config['assembly']['prefix'],vers),
        path=config['settings']['blobtools2_path'],
        covs=lambda wc: ' --cov '.join(["%s.%s.bam=%s" % (config['assembly']['prefix'], sra, sra) for sra in list_sra_accessions(reads)])
    conda:
        '../envs/blobtools2.yaml'
    threads: 1
    log:
      lambda wc: "logs/%s/blobtoolkit_replace_cov.log" % (config['assembly']['prefix'])
    resources:
        threads=1,
        btk=1
    shell:
        '{params.path}/blobtools replace \
            --cov {params.covs} \
            --threads {threads} \
            {params.id} > {log} 2>&1'


rule blobtoolkit_replace_busco:
    """
    import BUSCO results into BlobDir.
    """
    input:
        meta="%s%s/identifiers.json" % (config['assembly']['prefix'],vers),
        tsv=expand("%s_{lineage}.tsv" % config['assembly']['prefix'],lineage=config['busco']['lineages'])
    output:
        temp('busco.replaced'),
        expand("%s%s/{lineage}_busco.json" % (config['assembly']['prefix'],vers),lineage=config['busco']['lineages'])
    params:
        id="%s%s" % (config['assembly']['prefix'],vers),
        path=config['settings']['blobtools2_path'],
        busco=' --busco '.join(["%s_%s.tsv" % (config['assembly']['prefix'],lineage) for lineage in config['busco']['lineages']])
    conda:
        '../envs/blobtools2.yaml'
    threads: 1
    log:
      lambda wc: "logs/%s/blobtoolkit_replace_busco.log" % (config['assembly']['prefix'])
    resources:
        threads=1,
        btk=1
    shell:
        '{params.path}/blobtools replace \
            --busco {params.busco} \
            {params.id} > {log} 2>&1'
