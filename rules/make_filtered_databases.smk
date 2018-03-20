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
    conda:
         '../envs/py3.yaml'
    threads: 16
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
        '{name}.root.{root}{masked}/list'
    wildcard_constraints:
        root='\d+'
    params:
        mask_ids=lambda wc: similarity[wc.name]['mask_ids'],
        indir=lambda wc: "%s/split/%s" % (similarity[wc.name]['local'],wc.name)
    conda:
         '../envs/py3.yaml'
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
        lists='{name}.root.{root}{masked}/list'
    output:
        '{name}.root.{root}{masked}.dmnd'
    params:
        outfile=lambda wc: str("%s.root.%s%s.dmnd" % (wc.name,wc.root,wc.masked)),
        db=lambda wc: str("%s.root.%s%s" % (wc.name,wc.root,wc.masked)),
        tmpdir="%s" % config['settings']['tmp']
    wildcard_constraints:
        root='\d+'
    conda:
         '../envs/diamond.yaml'
    threads: 32
    resources:
        tmpdir=64
    priority:10
    shell:
        '{ENV} \
        mkdir -p {params.tmpdir} && \
        parallel --no-notice -m -P {threads} \
            "awk \'{{print \\$0}}\' {{}} | \
            awk \'{{print \\$0\\".fa\\"}}\' | \
            xargs cat" \
            :::: {input.lists} \
            > {params.tmpdir}/{params.db}.fa && \
        diamond makedb \
            --in {params.tmpdir}/{params.db}.fa \
            -p {threads} \
            -d {params.outfile} && \
        rm {params.tmpdir}/{params.db}.fa'


rule make_blast_db:
    """
    Generate a custom BLAST database from a list of per-taxon sequence
    files.
    """
    input:
        split=lambda wc: "%s/split/%s.done" % (similarity[wc.name]['local'],wc.name),
        idmap=lambda wc: "%s/full/%s.taxid_map.gz" % (similarity[wc.name]['local'],wc.name),
        lists='{name}.root.{root}{masked}/list'
    output:
        db='blast/{name}.root.{root}{masked}.{suffix}'
    wildcard_constraints:
        suffix='\wal',
        root='\d+'
    params:
        dbtype=lambda wc: 'prot' if wc.suffix == 'pal' else 'nucl',
        db=lambda wc: str("%s.root.%s%s" % (wc.name,wc.root,wc.masked)),
        tmpdir="%s" % config['settings']['tmp']
    conda:
         '../envs/blast.yaml'
    threads: 32
    resources:
        tmpdir=64
    shell:
        '{ENV} \
        set +o pipefail && \
        mkdir -p blast && \
        mkdir -p {params.tmpdir} && \
        > {params.tmpdir}/{params.db}.dblist && \
        parallel -j {threads} \
            --no-notice \
            "cat {{}} | \
            awk \'{{print \\$0\\".taxid_map\\"}}\' | \
            xargs cat > {params.tmpdir}/{params.db}.taxid_map_{{#}} && \
            cat {{}} | \
            awk \'{{print \\$0\\".fa\\"}}\' | \
            xargs cat > {params.tmpdir}/{params.db}.fa_{{#}} && \
            makeblastdb -dbtype {params.dbtype} \
                        -in {params.tmpdir}/{params.db}.fa_{{#}} \
                        -title {params.db}_{{#}} \
                        -out blast/{params.db}_{{#}} \
                        -parse_seqids \
                        -taxid_map {params.tmpdir}/{params.db}.taxid_map_{{#}} && \
            echo {params.db}_{{#}} >> {params.tmpdir}/{params.db}.dblist && \
            rm {params.tmpdir}/{params.db}.fa_{{#}} && \
            rm {params.tmpdir}/{params.db}.taxid_map_{{#}}" \
            :::: {input.lists} && \
        cd blast && \
        blastdb_aliastool -dblist_file {params.tmpdir}/{params.db}.dblist \
             -dbtype {params.dbtype} \
             -out {params.db} \
             -title {params.db} && \
        rm {params.tmpdir}/{params.db}.dblist'
