rule make_diamond_db:
    """
    Generate a custom Diamond database from a list of per-taxon sequence
    files.
    """
    input:
        split = lambda wc: "%s/split/%s.done" % (similarity[wc.name]['local'], wc.name),
        lists = 'blast/{name}.root.{root}{masked}.lists',
        nodes = "%s/nodes.dmp" % config['settings']['taxonomy']
    output:
        '{name}.root.{root}{masked}.dmnd'
    params:
        outfile = lambda wc: str("%s.root.%s%s.dmnd" % (wc.name, wc.root, wc.masked)),
        db = lambda wc: str("%s.root.%s%s" % (wc.name, wc.root, wc.masked)),
        indir = lambda wc: "%s/split/%s" % (similarity[wc.name]['local'], wc.name),
        tmpdir = "%s" % config['settings']['tmp']
    wildcard_constraints:
        root = r'\d+'
    conda:
        '../envs/diamond.yaml'
    threads: get_threads('make_diamond_db', maxcore)
    log:
        'logs/make_diamond_db/{name}.root.{root}{masked}.log'
    benchmark:
        'logs/make_diamond_db/{name}.root.{root}{masked}.benchmark.txt'
    resources:
        threads = get_threads('make_diamond_db', maxcore)
    shell:
        '(mkdir -p {params.tmpdir} && \
        echo "accession\taccession.version\ttaxid\tgi" > {params.tmpdir}/{params.db}.taxid_map && \
        parallel --no-notice -j {threads} \
            "gunzip -c {params.indir}/{{}}.taxid_map.gz" \
            :::: {input.lists} | \
                awk \'{{print $1 "\\t" $1 "\\t" $2 "\\t" 0}}\' >> \
                {params.tmpdir}/{params.db}.taxid_map && \
        parallel --no-notice -j {threads} \
            "seqtk subseq {params.indir}/{{}}.fa.gz blast/{params.db}_{{}}.accessions" \
            :::: {input.lists} | \
        diamond makedb \
            -p {threads} \
            -d {params.outfile} \
            --taxonmap {params.tmpdir}/{params.db}.taxid_map \
            --taxonnodes {input.nodes} && \
        rm {params.tmpdir}/{params.db}.taxid_map) 2> {log}'
