rule fetch_uniprot:
    """
    Fetch tarred Uniprot archive
    """
    output:
        temp("%s/full/{name}.tar.gz" % similarity['reference_proteomes']['local'])
    wildcard_constraints:
        # NB: the path to the local copy of the file must contain the string 'uniprot'
        path = '.+uniprot.+'
    params:
        ftp_url = 'ftp.ebi.ac.uk',
        ftp_dir = 'pub/databases/uniprot/current_release/knowledgebase/{name}'
    conda:
        '../envs/fetch.yaml'
    threads: get_threads('fetch_uniprot', 1)
    log:
        'logs/fetch_uniprot/{name}.log'
    benchmark:
        'logs/fetch_uniprot/{name}.benchmark.txt'
    resources:
        download = 1,
        threads = get_threads('fetch_uniprot', 1)
    shell:
        '(wget -q -O {output} \
            {params.ftp_url}/{params.ftp_dir}/$(curl \
            -vs {params.ftp_url}/{params.ftp_dir}/ 2>&1 | \
            awk \'/tar.gz/ {{print $9}}\')) 2> {log}'
