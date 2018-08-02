rule fetch_assembly:
    """
    Fetch a remote assembly from EBI or NCBI.
    """
    output:
        fa='{assembly}.fasta'
    params:
        url=lambda wc: prepare_ncbi_assembly_url(config['assembly']['accession'],config['assembly']['alias'])
    conda:
         '../envs/fetch.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'curl -s {params.url} | \
        pigz -d > {output.fa}'

rule fetch_sra:
    """
    Fetch sra archive from EBI.
    """
    output:
        temp('{sra}')
    params:
        sra = lambda wc: prepare_ebi_sra_url(wc.sra)
    wildcard_constraints:
        sra='\wRR\d+'
    conda:
         '../envs/fetch.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'wget -q {params.sra}'

rule extract_paired_fastq:
    """
    Extract paired fastq files from sra archive.
    """
    input:
        temp('{sra}')
    output:
        temp('{sra}_1.fastq.gz'),
        temp('{sra}_2.fastq.gz')
    conda:
         '../envs/fetch.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'fastq-dump --split-files {input} && \
        pigz {input}_*.fastq'

rule extract_fastq:
    """
    extract fastq files from sra archive.
    """
    input:
        temp('{sra}')
    output:
        temp('{sra}.fastq.gz')
    wildcard_constraints:
        sra='\wRR\d+'
    conda:
         '../envs/fetch.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'fastq-dump {input} && \
        pigz {input}.fastq'
