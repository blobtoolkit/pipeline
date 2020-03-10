rule split_fasta:
    """
    Use taxid_map to split FASTA into 100 files by taxid ending.
    Duplicate sequences in non-redundant databases so each taxon listed in
    header is represented by a single sequence.
    """
    input:
        fa = "%s/full/{name}.fa.gz" % similarity['reference_proteomes']['local'],
        idmap = "%s/full/{name}.taxid_map.gz" % similarity['reference_proteomes']['local']
    output:
        touch("%s/split/{name}.done" % similarity['reference_proteomes']['local'])
    params:
        tmpdir = lambda wc: "%s/%s" % (config['settings']['tmp'], wc.name),
        dir = uniprot_dir,
        chunk = config['settings']['chunk'],
        outdir = lambda wc: "%s/split/%s" % (similarity['reference_proteomes']['local'], wc.name),
    conda:
        '../envs/py3.yaml'
    threads: get_threads('split_fasta', multicore)
    log:
        'logs/expand_and_split_fasta/{name}.log'
    benchmark:
        'logs/expand_and_split_fasta/{name}.benchmark.txt'
    resources:
        tmpdir = 128
    script:
        '../scripts/expand_and_split_fasta.py'
