def list_similarity_results(config):
    """
    Generate a list of output filenames for sequence similarity searches
    based on list of databases in "config['similarity']".
    """
    path = []
    for db in config['similarity']['databases']:
        suffix = 'out' if db['tool'] == 'blast' else 'taxified.out'
        program = 'blastn' if db['type'] == 'nucl' else 'blastx' if db['tool'] == 'blast' else 'diamond'
        masked = ''
        if 'mask_ids' in db and isinstance(db['mask_ids'],(list,)):
            masked = "minus.%s" % '.'.join(str(mask) for mask in db['mask_ids'])
        else:
            masked = 'full'
        path.append("%s.%s.%s.root.%s.%s.%s" % (config['assembly']['name'],program,db['name'],db['root'],masked,suffix))
    return path

rule blobtools_create:
    """
    Use BlobTools create to generate a BlobDB using the ordered similarity
    searches and coverage files.
    """
    input:
        assembly='{assembly}.fna',
        dbs=list_similarity_results(config),
        coverage=expand('{{assembly}}.{sample}.bam.cov',sample=config['reads']['paired'])
    output:
        '{assembly}.blobDB.json'
    params:
        dbs=expand('-t {db}',db=list_similarity_results(config)),
        coverage=lambda wc: expand('-c '+wc.assembly+'.{sample}.bam.cov',sample=config['reads']['paired']),
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
