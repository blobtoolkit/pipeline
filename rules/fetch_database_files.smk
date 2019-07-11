rule fetch_ncbi_db:
    """
    Fetch formatted NCBI BLAST databases
    """
    output:
        '{path}/{name}.{suffix}'
    wildcard_constraints:
        # NB: the path to the local copy of the file must contain the string 'ncbi'
        path='.+ncbi.+',
        suffix='[np]al'
    params:
        ftp_url='ftp://ftp.ncbi.nlm.nih.gov',
        ftp_dir=lambda wc: "/blast/db%s" % '/v5' if wc.name.endswith('_v5') else ''
    conda:
         '../envs/fetch.yaml'
    threads: 1
    log:
      lambda wc: "logs/fetch_ncbi_db/%s.%s.log" % (wc.name, wc.suffix)
    resources:
        download=1,
        threads=1
    shell:
        '(wget "{params.ftp_url}{params.ftp_dir}{wildcards.name}.??.tar.gz" -P {wildcards.path}/ && \
        for file in {wildcards.path}/*.tar.gz; \
            do tar xf $file -C {wildcards.path} && rm $file; \
        done) 2> {log}'

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
        idmap=lambda wc: ncbi_idmap(wc.name)
    conda:
         '../envs/fetch.yaml'
    threads: 1
    log:
      lambda wc: "logs/fetch_ncbi_idmap/%s.log" % (wc.name)
    resources:
        download=1,
        threads=1
    shell:
        '(> {output.idmap} && \
        for x in "{params.idmap}"; do \
            wget -q -O idmap.gz $x \
            && zgrep -v taxid idmap.gz \
            | cut -f2,3 \
            | gzip >> {output.idmap} \
            && rm idmap.gz; \
        done) 2> {log}'

rule fetch_taxdump:
    """
    Fetch NCBI taxonomy "names.dmp" and "nodes.dmp" files
    """
    output:
        "%s/names.dmp" % config['settings']['taxonomy'],
        "%s/nodes.dmp" % config['settings']['taxonomy'],
        "%s/taxidlineage.dmp" % config['settings']['taxonomy']
    params:
        outdir=config['settings']['taxonomy']
    conda:
         '../envs/fetch.yaml'
    threads: 1
    log:
      "logs/fetch_taxdump.log"
    resources:
        download=1,
        threads=1
    shell:
        '(wget -q -O taxdump.tar.gz ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz \
        && tar xzvf taxdump.tar.gz \
        && mv *.dmp {params.outdir} \
        && rm taxdump.tar.gz) 2> {log}'

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
    log:
      lambda wc: "logs/fetch_uniprot/%s.log" % (wc.name)
    resources:
        download=1,
        threads=1
    shell:
        '(wget -q -O {output} \
            {params.ftp_url}/{params.ftp_dir}/$(curl \
            -vs {params.ftp_url}/{params.ftp_dir}/ 2>&1 | \
            awk \'/tar.gz/ {{print $9}}\')) 2> {log}'

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
    log:
      lambda wc: "logs/extract_uniprot/%s.log" % (wc.name)
    resources:
        tmpdir=24,
        threads=1
    script:
        '../scripts/extract_uniprot.py'

rule fetch_busco_lineage:
    """
    Fetch BUSCO lineages
    """
    output:
        gz=config['busco']['lineage_dir']+'/{lineage}.tar.gz',
        cfg=config['busco']['lineage_dir']+'/{lineage}/dataset.cfg'
    params:
        lineage=lambda wc: wc.lineage,
        dir=lambda wc: config['busco']['lineage_dir']
    conda:
         '../envs/fetch.yaml'
    threads: 1
    log:
      lambda wc: "logs/fetch_busco_lineage/%s.log" % (wc.lineage)
    resources:
        download=1,
        threads=1
    shell:
        '(wget -q -O {output.gz} "https://busco.ezlab.org/datasets/{params.lineage}.tar.gz" \
        && tar xf {output.gz} -C {params.dir}) 2> {log}'
