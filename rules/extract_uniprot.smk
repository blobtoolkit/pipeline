rule extract_uniprot:
    """
    Extract protein FASTA and idmapping files from Uniprot archive
    """
    input:
        "%s/full/{name}.tar.gz" % similarity['reference_proteomes']['local']
    output:
        "%s/full/{name}.fa.gz" % similarity['reference_proteomes']['local'],
        "%s/full/{name}.taxid_map.gz" % similarity['reference_proteomes']['local']
    params:
        tmpdir = lambda wc: "%s/%s" % (config['settings']['tmp'], wc.name),
        dir = uniprot_dir,
    conda:
        '../envs/py3.yaml'
    threads: get_threads('extract_uniprot', 1)
    log:
        'logs/extract_uniprot/{name}.log'
    benchmark:
        'logs/extract_uniprot/{name}.benchmark.txt'
    resources:
        tmpdir = 24,
        threads = get_threads('extract_uniprot', 1)
    script:
        '../scripts/extract_uniprot.py'
