rule fetch_ncbi_fasta:
    """
    Fetch FASTA files corresponding to NCBI BLAST databases
    """
    output:
        fa='{path}/full/{name}.fa.gz'
    wildcard_constraints:
        # NB: the path to the local copy of the file must contain the string 'ncbi'
        path='.+ncbi.+'
    params:
        ftp_url='ftp.ncbi.nlm.nih.gov',
        ftp_dir='/blast/db/FASTA'
    conda:
         '../envs/fetch.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'curl {params.ftp_url}{params.ftp_dir}/{wildcards.name}.gz > {output.fa}'

rule fetch_ncbi_idmap:
    """
    Fetch FASTA files corresponding to NCBI BLAST databases
    """
    output:
        idmap='{path}/full/{name}.taxid_map.gz'
    wildcard_constraints:
        # NB: the path to the local copy of the file must contain the string 'ncbi'
        path='.+ncbi.+'
    params:
        idmap=lambda wildcards: ncbi_idmap(wildcards.name)
    conda:
         '../envs/fetch.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        '> {output.idmap} && \
        for x in "{params.idmap}"; do \
            curl $x | zgrep -v taxid | cut -f2,3 | gzip >> {output.idmap}; \
        done'

rule fetch_taxdump:
    """
    Fetch NCBI taxonomy "names.dmp" and "nodes.dmp" files
    """
    output:
        "%s/names.dmp" % config['settings']['taxonomy'],
        "%s/nodes.dmp" % config['settings']['taxonomy'],
        temp('citations.dmp'),
        temp('delnodes.dmp'),
        temp('division.dmp'),
        temp('gencode.dmp'),
        temp('merged.dmp'),
        temp('gc.prt')
    params:
        outdir=config['settings']['taxonomy']
    conda:
         '../envs/fetch.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'curl -vs ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz \
        | tar xzvf - \
        && mv n*s.dmp {params.outdir}'

# rule update_blobtools_nodesdb:
#     """
#     Update the BlobTools "nodesDB.txt" file
#     """
#     input:
#         names="%s/names.dmp" % config['settings']['taxonomy'],
#         nodes="%s/nodes.dmp" % config['settings']['taxonomy']
#     output:
#         "%s/data/nodesDB.txt" % config['settings']['blobtools']
#     conda:
#          '../envs/blobtools.yaml'
#     threads: 1
#     shell:
#         '{ENV} \
#         PATH=' + config['settings']['blobtools'] + ':$PATH && \
#         blobtools nodesdb \
#             --nodes {input.nodes} \
#             --names {input.names}'

rule fetch_uniprot:
    """
    Fetch tarred Uniprot archive
    """
    output:
        temp('{path}/full/{name}.tar.gz')
    wildcard_constraints:
        # NB: the path to the local copy of the file must contain the string 'uniprot'
        path='.+uniprot.+'
    params:
        ftp_url='ftp.ebi.ac.uk',
        ftp_dir='pub/databases/uniprot/current_release/knowledgebase/{name}'
    conda:
         '../envs/fetch.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'curl {params.ftp_url}/{params.ftp_dir}/$(curl \
            -vs {params.ftp_url}/{params.ftp_dir}/ 2>&1 | \
            awk \'/tar.gz/ {{print $9}}\') > {output}'

rule extract_uniprot:
    """
    Extract protein FASTA and idmapping files from Uniprot archive
    """
    input:
        '{path}/full/{name}.tar.gz'
    output:
        '{path}/full/{name}.fa.gz',
        '{path}/full/{name}.taxid_map.gz'
    params:
        tmpdir=lambda wc: "%s/%s" % (config['settings']['tmp'],wc.name)
    wildcard_constraints:
        # NB: the path to the local copy of the file must contain the string 'uniprot'
        path='.+uniprot.+'
    conda:
         '../envs/py3.yaml'
    threads: 1
    resources:
        tmpdir=24,
        threads=1
    script:
        '../scripts/extract_uniprot.py'
