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
    benchmark:
        'logs/{assembly}/generate_metadata.benchmark.txt'
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
        "{assembly}%s/identifiers.json" % rev
    params:
        taxdump=config['settings']['taxonomy'],
        id=lambda wc: "%s%s" % (wc.assembly,rev),
        path=config['settings']['blobtools2_path'],
        taxid=config['taxon']['taxid'],
    conda:
        '../envs/blobtools2.yaml'
    threads: 1
    log:
        lambda wc: "logs/%s/blobtoolkit_create.log" % (wc.assembly)
    benchmark:
        'logs/{assembly}/blobtoolkit_create.benchmark.txt'
    resources:
        threads=1,
        btk=1
    shell:
        '{params.path}/blobtools replace \
            --fasta {input.fasta} \
            --meta {input.yaml} \
            --taxdump "{params.taxdump}" \
            --taxid {params.taxid} \
            {params.id} > {log} 2>&1'

rule blobtoolkit_add_hits:
    """
    Add ordered similarity search results to a BlobDir.
    """
    input:
        meta="%s%s/identifiers.json" % (config['assembly']['prefix'],rev),
        dbs=list_similarity_results(config),
        lineages="%s/taxidlineage.dmp" % (config['settings']['taxonomy'])
    output:
        "{assembly}%s/%s_phylum_positions.json" % (rev,config['similarity']['taxrule'])
    params:
        taxrule=config['similarity']['taxrule'] if 'taxrule' in config['similarity'] else 'bestsumorder',
        taxdump=config['settings']['taxonomy'],
        id=lambda wc: "%s%s" % (wc.assembly,rev),
        path=config['settings']['blobtools2_path'],
        dbs='.raw --hits '.join(list_similarity_results(config))
    conda:
        '../envs/blobtools2.yaml'
    threads: 1
    log:
        lambda wc: "logs/%s/blobtoolkit_add_hits.log" % (wc.assembly)
    benchmark:
        'logs/{assembly}/blobtoolkit_add_hits.benchmark.txt'
    resources:
        threads=1,
        btk=1
    shell:
        '{params.path}/blobtools replace \
            --hits {params.dbs} \
            --taxrule "{params.taxrule}" \
            --taxdump "{params.taxdump}" \
            {params.id} > {log} 2>&1'

rule blobtoolkit_add_cov:
    """
    Use BlobTools2 add to add coverage to a BlobDir from BAM files.
    """
    input:
        meta="%s%s/identifiers.json" % (config['assembly']['prefix'],rev),
        bam=expand("%s.{sra}.bam" % asm, sra=list_sra_accessions(reads))
    output:
        expand("%s%s/{sra}_cov.json" % (config['assembly']['prefix'],rev),sra=list_sra_accessions(reads))
    params:
        id="%s%s" % (config['assembly']['prefix'],rev),
        path=config['settings']['blobtools2_path'],
        covs=lambda wc: ' --cov '.join(["%s.%s.bam=%s" % (config['assembly']['prefix'], sra, sra) for sra in list_sra_accessions(reads)])
    conda:
        '../envs/blobtools2.yaml'
    threads: 1
    log:
        lambda wc: "logs/%s/blobtoolkit_add_cov.log" % (config['assembly']['prefix'])
    benchmark:
        "logs/%s/blobtoolkit_add_cov.benchmark.txt" % (config['assembly']['prefix'])
    resources:
        threads=1,
        btk=1
    shell:
        '{params.path}/blobtools replace \
            --cov {params.covs} \
            --threads {threads} \
            {params.id} > {log} 2>&1'


rule blobtoolkit_add_busco:
    """
    import BUSCO results into BlobDir.
    """
    input:
        meta="%s%s/identifiers.json" % (config['assembly']['prefix'],rev),
        tsv=expand("%s_{lineage}.tsv" % config['assembly']['prefix'],lineage=config['busco']['lineages'])
    output:
        expand("%s%s/{lineage}_busco.json" % (config['assembly']['prefix'],rev),lineage=config['busco']['lineages'])
    params:
        id="%s%s" % (config['assembly']['prefix'],rev),
        path=config['settings']['blobtools2_path'],
        busco=' --busco '.join(["%s_%s.tsv" % (config['assembly']['prefix'],lineage) for lineage in config['busco']['lineages']])
    conda:
        '../envs/blobtools2.yaml'
    threads: 1
    log:
        lambda wc: "logs/%s/blobtoolkit_add_busco.log" % (config['assembly']['prefix'])
    benchmark:
        "logs/%s/blobtoolkit_add_busco.benchmark.txt" % (config['assembly']['prefix'])
    resources:
        threads=1,
        btk=1
    shell:
        '{params.path}/blobtools replace \
            --busco {params.busco} \
            {params.id} > {log} 2>&1'
