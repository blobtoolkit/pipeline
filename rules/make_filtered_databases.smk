rule expand_and_split_fasta:
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
    threads: 8
    resources:
        tmpdir=128,
        threads=8
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
        '{name}.root.{root}{masked}.taxids'
    wildcard_constraints:
        root='\d+'
    params:
        mask_ids=lambda wc: similarity[wc.name]['mask_ids'],
        db=lambda wc: str("%s.root.%s%s" % (wc.name,wc.root,wc.masked))
    conda:
         '../envs/py3.yaml'
    threads: 1
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
        nodes="%s/nodes.dmp" % config['settings']['taxonomy'],
        split=lambda wc: "%s/split/%s.done" % (similarity[wc.name]['local'],wc.name)
    output:
        'blast/{name}.root.{root}{masked}.lists'
    wildcard_constraints:
        root='\d+'
    params:
        mask_ids=lambda wc: similarity[wc.name]['mask_ids'],
        db=lambda wc: str("%s.root.%s%s" % (wc.name,wc.root,wc.masked)),
        indir=lambda wc: "%s/split/%s" % (similarity[wc.name]['local'],wc.name),
        chunk=config['settings']['chunk']
    conda:
         '../envs/py3.yaml'
    threads: 1
    resources:
        threads=1
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
    threads: 32
    resources:
        threads=32
    shell:
        'mkdir -p {params.tmpdir} && \
        parallel --no-notice -j {threads} \
            "gunzip {params.indir}/{{}}.taxid_map.gz" \
            :::: {input.lists} > \
                {params.tmpdir}/{params.db}.taxid_map && \
        parallel --no-notice -j {threads} \
            "seqtk subseq {params.indir}/{{}}.fa.gz blast/{params.db}_{{}}.accessions" \
            :::: {input.lists} | \
        diamond makedb \
            -p {threads} \
            -d {params.outfile} \
            --taxonmap {params.tmpdir}/{params.db}.taxid_map \
            --taxonnodes {input.nodes} && \
        rm {params.tmpdir}/{params.db}.taxid_map'


# rule make_blast_db:
#     """
#     Generate a custom BLAST database from a list of per-taxon sequence
#     files.
#     """
#     input:
#         split=lambda wc: "%s/split/%s.done" % (similarity[wc.name]['local'],wc.name),
#         lists='blast/{name}.root.{root}{masked}.lists'
#     output:
#         db='blast/{name}.root.{root}{masked}.{suffix}'
#     wildcard_constraints:
#         suffix='\wal',
#         root='\d+'
#     params:
#         dbtype=lambda wc: 'prot' if wc.suffix == 'pal' else 'nucl',
#         indir=lambda wc: "%s/split/%s" % (similarity[wc.name]['local'],wc.name),
#         db=lambda wc: str("%s.root.%s%s" % (wc.name,wc.root,wc.masked)),
#         tmpdir="%s" % config['settings']['tmp']
#     conda:
#          '../envs/blast.yaml'
#     threads: 32
#     resources:
#         tmpdir=64,
#         threads=32
#     shell:
#         'set +o pipefail && \
#         mkdir -p blast && \
#         mkdir -p {params.tmpdir} && \
#         > {params.tmpdir}/{params.db}.dblist && \
#         parallel -j {threads} \
#             --no-notice \
#             "pigz -dc {params.indir}/{{}}.taxid_map.gz \
#                 > {params.tmpdir}/{params.db}_{{}}.taxid_map && \
#             seqtk subseq {params.indir}/{{}}.fa.gz blast/{params.db}_{{}}.accessions | \
#             makeblastdb -dbtype {params.dbtype} \
#                         -title {params.db}_{{}} \
#                         -out blast/{params.db}_{{}} \
#                         -parse_seqids \
#                         -taxid_map {params.tmpdir}/{params.db}_{{}}.taxid_map && \
#             echo {params.db}_{{}} >> {params.tmpdir}/{params.db}.dblist && \
#             rm {params.tmpdir}/{params.db}_{{}}.taxid_map" \
#             :::: {input.lists} && \
#         cd blast && \
#         blastdb_aliastool -dblist_file {params.tmpdir}/{params.db}.dblist \
#              -dbtype {params.dbtype} \
#              -out {params.db} \
#              -title {params.db} && \
#         rm {params.tmpdir}/{params.db}.dblist'
