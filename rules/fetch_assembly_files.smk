rule fetch_assembly:
    """
    Fetch a remote assembly from EBI or NCBI.
    """
    output:
        fa='{assembly}.fasta'
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
        temp('{sra}_1.fastq.gz'),
        temp('{sra}_2.fastq.gz')
    params:
        sra = lambda wc: wc.sra
    conda:
         '../envs/py3.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'enaDataGet -f fastq {params.sra} && \
        mv {params.sra}/{params.sra}_1.fastq.gz {output[0]} && \
        mv {params.sra}/{params.sra}_2.fastq.gz {output[1]} && \
        rm -r {params.sra}'

rule fetch_fastq:
    """
    Fetch paired fastq files from EBI.
    """
    output:
        temp('{sra}.fastq.gz')
    params:
        sra = lambda wc: wc.sra
    wildcard_constraints:
        sra='\wRR\d+'
    conda:
         '../envs/py3.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'enaDataGet -f fastq {params.sra} && \
        rename -f "s/(_subreads|_consensus)//" {params.sra}/*.fastq.gz && \
        mv {params.sra}/{params.sra}.fastq.gz {output[1]} && \
        rm -r {params.sra}'
