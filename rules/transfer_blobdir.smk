#    input:
#        "{asm}.done"

rule validate_dataset:
    """
    Run BlobToolKit validator on a dataset to check all expected fields are present.
    """
    input:
        cov=expand("%s/{sra}_cov.json" % asm,sra=list_sra_accessions(reads)),
        tax="%s/%s_phylum_positions.json" % (asm,config['similarity']['taxrule']),
        busco=expand("%s/{lineage}_busco.json" % asm,lineage=config['busco']['lineages']),
        ids="%s/identifiers.json" % asm
    output:
        '{assembly}.valid'
    params:
        assembly = lambda wc: wc.assembly
    conda:
         '../envs/blobtools2.yaml'
    threads: 1
    resources:
        threads=1
    shell:
        'pip install fastjsonschema \
        && /ceph/software/blobtoolkit/specification/validate.py {params.assembly}/meta.json \
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
    Use BTK CLI to generate a set of static images.
    """
    input:
        valid='{assembly}.valid',
        hosted='{assembly}.hosted',
        cov=expand("%s/{sra}_cov.json" % asm,sra=list_sra_accessions(reads))
    output:
        '{assembly}/summary.json'
    params:
        assembly=lambda wc: wc.assembly,
        port=8080
    conda:
         '../envs/blobtools2.yaml'
    threads: 1
    resources:
        threads=1
    script:
        '../scripts/generate_static_images.py'
