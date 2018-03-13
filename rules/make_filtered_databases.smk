# ===================================================================== #
# NB: The rule split_fasta_by_taxid takes too long to run.                  #
# Need to improve efficiency but size of dataset caused problems with a #
# simple multiprocessing approach.                                      #
# ===================================================================== #

rule split_fasta_by_taxid:
    """
    Use taxid_map to split FASTA into one file per taxon.
    Duplicate sequences in non-redundant databases so each taxon listed in
    header is represented by a single sequence.
    """
    input:
        fa='{path}/full/{name}.fa.gz',
        idmap='{path}/full/{name}.taxid_map.gz'
    output:
        touch('{path}/split/{name}.done')
    params:
        tmpdir="%s/by_taxid" % config['settings']['tmp'],
        chunk=config['settings']['chunk']
    threads: 1
    resources:
        tmpdir=128
    script:
        '../scripts/split_fasta_by_taxid.py'

rule list_sequence_files:
    """
    Generate a list of per-taxon sequence files needed to create a custom
    database containing all descendants of a specified root, optionally
    with one or more lineages masked.
    """
    input:
        nodes="%s/nodes.dmp" % config['settings']['taxonomy'],
        split=lambda wc: "%s/split/%s.done" % (similarity[wc.name]['local'],wc.name)
    output:
        '{name}.root.{root}{masked}list.gz'
    wildcard_constraints:
        root='\d+'
    params:
        mask_ids=lambda wc: similarity[wc.name]['mask_ids'],
        indir=lambda wc: "%s/split/%s" % (similarity[wc.name]['local'],wc.name)
    threads: 1
    script:
        '../scripts/masked_list_by_root.py'

rule make_diamond_db:
    """
    Generate a custom Diamond database from a list of per-taxon sequence
    files.
    """
    input:
        split=lambda wc: "%s/split/%s.done" % (similarity[wc.name]['local'],wc.name),
        lists='{name}.root.{root}{masked}list.gz'
    output:
        '{name}.root.{root}{masked}dmnd'
    params:
        outfile=lambda wc: str("%s.root.%s%sdmnd" % (wc.name,wc.root,wc.masked)).rstrip('.')
    wildcard_constraints:
        root='\d+'
    threads: 32
    resources:
        tmpdir=64
    priority:10
    shell:
        '{ENV} \
        pigz -dc {input.lists} | \
        xargs cat | \
        pigz -dc | \
        diamond makedb \
            -p {threads} \
            -d {params.outfile}'


rule make_blast_db:
    """
    Generate a custom BLAST database from a list of per-taxon sequence
    files.
    """
    input:
        split=lambda wc: "%s/split/%s.done" % (similarity[wc.name]['local'],wc.name),
        idmap=lambda wc: "%s/full/%s.taxid_map.gz" % (similarity[wc.name]['local'],wc.name),
        lists='{name}.root.{root}{masked}list.gz'
    output:
        db='{name}.root.{root}{masked}{suffix}'
    wildcard_constraints:
        suffix='\wal',
        root='\d+'
    params:
        dbtype=lambda wc: 'prot' if wc.suffix == 'pal' else 'nucl',
        db=lambda wc: str("%s.root.%s%s" % (wc.name,wc.root,wc.masked)).rstrip('.'),
        block='500M'
    threads: 32
    resources:
        tmpdir=64
    shell:
        '{ENV} \
        set +o pipefail && \
        > {params.db}.dblist && \
        pigz -dc {input.lists} | \
        xargs cat | \
        pigz -dc | \
        parallel -k \
            -j {threads} \
            --block {params.block} \
            --recstart \'>\' \
            --no-notice \
            --pipe \
            "makeblastdb -dbtype {params.dbtype} \
                         -title {params.db}_{{#}} \
                         -out {params.db}_{{#}} \
                         -parse_seqids \
                         -taxid_map <(pigz -dc {input.idmap}) && \
            echo {params.db}_{{#}} >> {params.db}.dblist" && \
        blastdb_aliastool -dblist_file {params.db}.dblist \
             -dbtype {params.dbtype} \
             -out {params.db} \
             -title {params.db} && \
        rm {params.db}.dblist'
