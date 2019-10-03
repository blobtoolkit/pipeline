rule validate_dataset:
    """
    Run BlobToolKit validator on a dataset to check all expected fields are present.
    """
    input:
        cov=expand("%s%s/{sra}_cov.json" % (asm,rev),sra=list_sra_accessions(reads)),
        tax="%s%s/%s_phylum_positions.json" % (asm,rev,config['similarity']['taxrule']),
        busco=expand("%s%s/{lineage}_busco.json" % (asm,rev),lineage=config['busco']['lineages']),
        ids="%s%s/identifiers.json" % (asm,rev)
    output:
        temp('{assembly}.valid')
    params:
        assembly = lambda wc: wc.assembly
    conda:
         '../envs/blobtools2.yaml'
    threads: 1
    log:
      lambda wc: "logs/%s/validate_dataset.log" % (wc.assembly)
    resources:
        threads=1
    shell:
        'validate.py {params.assembly}/meta.json > {log} 2>&1 \
        && touch {params.assembly}.valid'


rule host_dataset:
    """
    Host dataset with BTK Viewer.
    """
    input:
        valid='{assembly}.valid'
    output:
        temp('{assembly}.hosted')
    params:
        assembly=lambda wc: wc.assembly,
        port=8080
    conda:
         '../envs/blobtools2.yaml'
    threads: 1
    resources:
        threads=1
    shell:
        'PID="$(lsof -nP -iTCP:{params.port} | grep LISTEN | awk \'{{print $2}}\')"; \
        PID=($PID); \
        [[ "${{PID:-}}" == "-" ]] && {{ echo "BlobToolKit Viewer is not running on port {params.port}"; exit 1; }}; \
        touch {params.assembly}.hosted'


rule generate_images:
    """
    Use BlobTools2 view to generate a set of static images.
    """
    input:
        valid='{assembly}.valid',
        hosted='{assembly}.hosted',
        cov=expand("%s%s/{sra}_cov.json" % (asm,rev),sra=list_sra_accessions(reads))
    output:
        '{assembly}/cumulative.png'
    params:
        assembly=lambda wc: wc.assembly,
        host='http://localhost:8080'
    conda:
         '../envs/blobtools2.yaml'
    threads: 1
    log:
      lambda wc: "logs/%s/generate_images.log" % (wc.assembly)
    resources:
        threads=1
    script:
        '../scripts/generate_static_images.py'


rule generate_summary:
    """
    Use BlobTools2 filter to generate a dataset summary.
    """
    input:
        '{assembly}/cumulative.png'
    output:
        '{assembly}/summary.json'
    params:
        assembly=lambda wc: wc.assembly
    conda:
         '../envs/blobtools2.yaml'
    threads: 1
    log:
      lambda wc: "logs/%s/generate_summary.log" % (wc.assembly)
    resources:
        threads=1
    shell:
        'blobtools filter --summary {output} {params.assembly} 2> {log}'


rule checksum_files:
    """
    Calculate SHA1 checksum for all files in dataset.
    """
    input:
        '{assembly}/summary.json'
    output:
        '{assembly}/CHECKSUM'
    params:
        assembly=lambda wc: wc.assembly
    threads: 1
    log:
      lambda wc: "logs/%s/checksum_files.log" % (wc.assembly)
    resources:
        threads=1
    shell:
        '(find {params.assembly}/ -type f -exec sha1sum {{}} \';\' \
        | sort -k 2 \
        | sed \'s/{params.assembly}\///\' > {params.assembly}/CHECKSUM) 2> {log}'


rule transfer_dataset:
    """
    Transfer dataset out of working directory.
    """
    input:
        '{assembly}/CHECKSUM'
    output:
        temp('{assembly}.complete'),
    params:
        assembly=lambda wc: wc.assembly,
        destdir=config['destdir'],
    threads: 1
    log:
      lambda wc: "logs/%s/transfer_dataset.log" % (wc.assembly)
    resources:
        threads=1
    shell:
        '(rsync -av --remove-source-files {params.assembly}* {params.destdir}/ \
        && rmdir {params.assembly} \
        && touch {params.assembly}.complete) 2> {log}'
