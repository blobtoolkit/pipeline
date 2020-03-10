def hits_inputs(config, rev):
    return {
        'meta': "%s%s/identifiers.json" % (config['assembly']['prefix'], rev),
        'dbs': list_similarity_results(config),
        'lineages': "%s/taxidlineage.dmp" % (config['settings']['taxonomy'])
    }


def hits_params(config, rev, aa=False):
    taxrule = config['similarity']['taxrule']
    files = list_similarity_results(config)
    if aa:
        files.reverse()
    return {
        'taxrule': taxrule,
        'taxdump': taxdump_dir,
        'id': lambda wc: "%s%s" % (wc.assembly, rev),
        'dbs': '.raw --hits '.join(files),
        'path': config['settings']['blobtools2_path']
    }


rule blobtoolkit_add_hits:
    """
    Add ordered similarity search results to a BlobDir.
    """
    input:
        **hits_inputs(config, rev)
    output:
        "{assembly}%s/%s_phylum_positions.json" % (rev, config['similarity']['taxrule'])
    params:
        **hits_params(config, rev)
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('blobtoolkit_add_hits', 1)
    log:
        'logs/{assembly}/blobtoolkit_add_hits.log'
    benchmark:
        'logs/{assembly}/blobtoolkit_add_hits.benchmark.txt'
    resources:
        threads = get_threads('blobtoolkit_add_hits', 1),
        btk = 1
    shell:
        'PATH={params.path}:$PATH && \
        blobtools replace \
            --hits {params.dbs} \
            --taxrule "{params.taxrule}" \
            --taxdump "{params.taxdump}" \
            {params.id} > {log} 2>&1'


rule blobtoolkit_add_nt_hits:
    """
    Add ordered similarity search results to a BlobDir.
    """
    input:
        **hits_inputs(config, rev)
    output:
        "{assembly}%s/%s_nt_phylum_positions.json" % (rev, config['similarity']['taxrule'].replace('each', 'best'))
    params:
        **hits_params(config, rev)
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('blobtoolkit_add_hits', 1)
    log:
        'logs/{assembly}/blobtoolkit_add_nt_hits.log'
    benchmark:
        'logs/{assembly}/blobtoolkit_add_nt_hits.benchmark.txt'
    resources:
        threads = get_threads('blobtoolkit_add_hits', 1),
        btk = 1
    shell:
        'PATH={params.path}:$PATH && \
        blobtools replace \
            --hits {params.dbs} \
            --taxrule "{params.taxrule}={params.taxrule}_nt" \
            --taxdump "{params.taxdump}" \
            {params.id} > {log} 2>&1'


rule blobtoolkit_add_aa_hits:
    """
    Add ordered similarity search results to a BlobDir.
    """
    input:
        **hits_inputs(config, rev)
    output:
        "{assembly}%s/%s_aa_phylum_positions.json" % (rev, config['similarity']['taxrule'].replace('each', 'best'))
    params:
        **hits_params(config, rev, True)
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('blobtoolkit_add_hits', 1)
    log:
        'logs/{assembly}/blobtoolkit_add_aa_hits.log'
    benchmark:
        'logs/{assembly}/blobtoolkit_add_aa_hits.benchmark.txt'
    resources:
        threads = get_threads('blobtoolkit_add_hits', 1),
        btk = 1
    priority:
        100
    shell:
        'PATH={params.path}:$PATH && \
        blobtools replace \
            --hits {params.dbs} \
            --taxrule "{params.taxrule}={params.taxrule}_aa" \
            --taxdump "{params.taxdump}" \
            {params.id} > {log} 2>&1'
