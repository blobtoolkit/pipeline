rule make_uniprot_db:
    """
    Generate a uniprot diamond database.
    """
    input:
        fasta = "%s/full/{name}.fa.gz" % similarity['reference_proteomes']['local'],
        idmap = "%s/full/{name}.taxid_map.gz" % similarity['reference_proteomes']['local'],
        nodes = "%s/nodes.dmp" % config['settings']['taxonomy']
    output:
        "%s/full/{name}.dmnd" % similarity['reference_proteomes']['local']
    params:
        db = lambda wc: wc.name,
        tmpdir = "%s" % config['settings']['tmp']
    conda:
        '../envs/diamond.yaml'
    threads: get_threads('make_uniprot_db', maxcore)
    log:
        'logs/{name}.make_uniprot_db.log'
    benchmark:
        'logs/{name}.make_uniprot_db.benchmark.txt'
    resources:
        threads = get_threads('make_uniprot_db', maxcore)
    shell:
        '(mkdir -p {params.tmpdir} && \
        echo "accession\taccession.version\ttaxid\tgi" > {params.tmpdir}/{params.db}.taxid_map && \
        gunzip -c {input.idmap} | \
                awk \'{{print $1 "\\t" $1 "\\t" $2 "\\t" 0}}\' > \
                {params.tmpdir}/{params.db}.taxid_map && \
        gunzip -c {input.fasta} | \
        diamond makedb \
            -p {threads} \
            -d {output} \
            --taxonmap {params.tmpdir}/{params.db}.taxid_map \
            --taxonnodes {input.nodes} && \
        rm {params.tmpdir}/{params.db}.taxid_map) 2> {log}'
