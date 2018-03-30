rule fetch_assembly:
    """
    Fetch a remote assembly from EBI or NCBI.
    """
    output:
        fa=temp('{assembly}.fasta')
    params:
        assembly=lambda wc: wc.assembly
    conda:
         '../envs/py3.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'enaDataGet -f fasta {params.assembly} && \
        pigz -d {output}.gz'

rule fetch_paired_fastq:
    """
    Fetch paired fastq files from EBI.
    """
    output:
        temp('{sra}/{sra}_1.fastq.gz'),
        temp('{sra}/{sra}_2.fastq.gz')
    params:
        sra = lambda wc: wc.sra
    conda:
         '../envs/py3.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'enaDataGet -f fastq {params.sra}'

rule fetch_fastq:
    """
    Fetch paired fastq files from EBI.
    """
    output:
        temp('{sra}/{sra}.fastq.gz')
    params:
        sra = lambda wc: wc.sra
    wildcard_constraints:
        sra='\w\w\w\d+'
    conda:
         '../envs/py3.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'enaDataGet -f fastq {params.sra}'
