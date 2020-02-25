rule blobtoolkit_add_hits:
    """
    Add ordered similarity search results to a BlobDir.
    """
    input:
        meta = "%s%s/identifiers.json" % (config['assembly']['prefix'], rev),
        dbs = list_similarity_results(config),
        lineages = "%s/taxidlineage.dmp" % (config['settings']['taxonomy'])
    output:
        "{assembly}%s/%s_phylum_positions.json" % (rev, config['similarity']['taxrule'])
    params:
        taxrule = config['similarity']['taxrule'] if 'taxrule' in config['similarity'] else 'bestsumorder',
        taxdump = config['settings']['taxonomy'],
        id = lambda wc: "%s%s" % (wc.assembly, rev),
        path = config['settings']['blobtools2_path'],
        dbs = '.raw --hits '.join(list_similarity_results(config))
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
        '{params.path}/blobtools replace \
            --hits {params.dbs} \
            --taxrule "{params.taxrule}" \
            --taxdump "{params.taxdump}" \
            {params.id} > {log} 2>&1'
