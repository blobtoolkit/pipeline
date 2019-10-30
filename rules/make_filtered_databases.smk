rule split_fasta:
    """
    Use taxid_map to split FASTA into 100 files by taxid ending.
    Duplicate sequences in non-redundant databases so each taxon listed in
    header is represented by a single sequence.
    """
    input:
        fa='{path}/full/{name}.fa.gz',
        idmap='{path}/full/{name}.taxid_map.gz'
    output:
        touch('{path}/split/{name}.done')
    params:
        tmpdir=lambda wc: "%s/%s" % (config['settings']['tmp'],wc.name),
        chunk=config['settings']['chunk'],
        outdir=lambda wc: "%s/split/%s" % (wc.path,wc.name),
    conda:
        '../envs/py3.yaml'
    threads: lambda x: multicore
    log:
        lambda wc: "logs/expand_and_split_fasta/%s.log" % (wc.name)
    benchmark:
        'logs/expand_and_split_fasta/{name}.benchmark.txt'
    resources:
        tmpdir=128,
        threads=lambda x: multicore
    script:
        '../scripts/expand_and_split_fasta.py'

rule make_taxid_list:
    """
    Generate a list of taxids containing all descendants of a specified root,
    optionally with one or more lineages masked.
    """
    input:
        nodes="%s/nodes.dmp" % config['settings']['taxonomy']
    output:
        '{name}.root.{root}{masked}.taxids',
        '{name}.root.{root}{masked}.negative.taxids'
    wildcard_constraints:
        root='\d+'
    params:
        mask_ids=lambda wc: similarity[wc.name]['mask_ids'],
        db=lambda wc: str("%s.root.%s%s" % (wc.name,wc.root,wc.masked))
    conda:
        '../envs/py3.yaml'
    threads: 1
    log:
        lambda wc: "logs/make_taxid_list/%s.root.%s%s.log" % (wc.name, wc.root, wc.masked)
    benchmark:
        'logs/make_taxid_list/{name}.root.{root}{masked}.benchmark.txt'
    resources:
        threads=1
    script:
        '../scripts/make_taxid_list.py'

rule make_masked_lists:
    """
    Generate a list of accessions needed to create a custom
    database containing all descendants of a specified root, optionally
    with one or more lineages masked.
    """
    input:
        split=lambda wc: "%s/split/%s.done" % (similarity[wc.name]['local'],wc.name),
        taxids='{name}.root.{root}{masked}.negative.taxids'
    output:
        'blast/{name}.root.{root}{masked}.lists'
    wildcard_constraints:
        root='\d+'
    params:
        db=lambda wc: str("%s.root.%s%s" % (wc.name,wc.root,wc.masked)),
        indir=lambda wc: "%s/split/%s" % (similarity[wc.name]['local'],wc.name),
        chunk=config['settings']['chunk']
    conda:
        '../envs/py3.yaml'
    threads: lambda x: maxcore
    log:
        lambda wc: "logs/make_masked_lists/%s.root.%s%s.log" % (wc.name, wc.root, wc.masked)
    benchmark:
        lambda wc: "logs/make_masked_lists/{name}.root.{root}{masked}.benchmark.txt'
    resources:
        threads=lambda x: maxcore
    script:
        '../scripts/make_masked_lists.py'

rule make_diamond_db:
    """
    Generate a custom Diamond database from a list of per-taxon sequence
    files.
    """
    input:
        split=lambda wc: "%s/split/%s.done" % (similarity[wc.name]['local'],wc.name),
        lists='blast/{name}.root.{root}{masked}.lists',
        nodes="%s/nodes.dmp" % config['settings']['taxonomy']
    output:
        '{name}.root.{root}{masked}.dmnd'
    params:
        outfile=lambda wc: str("%s.root.%s%s.dmnd" % (wc.name,wc.root,wc.masked)),
        db=lambda wc: str("%s.root.%s%s" % (wc.name,wc.root,wc.masked)),
        indir=lambda wc: "%s/split/%s" % (similarity[wc.name]['local'],wc.name),
        tmpdir="%s" % config['settings']['tmp']
    wildcard_constraints:
        root='\d+'
    conda:
        '../envs/diamond.yaml'
    threads: lambda x: maxcore
    log:
        lambda wc: "logs/make_diamond_db/%s.root.%s%s.log" % (wc.name, wc.root, wc.masked)
    benchmark:
        'logs/make_diamond_db/{name}.root.{root}{masked}.benchmark.txt'
    resources:
        threads=lambda x: maxcore
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
