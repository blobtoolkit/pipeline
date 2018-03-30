# def assembly_url(prefix):
#     """
#     Generate a URL for an assembly at EBI.
#     """
#     return "https://www.ebi.ac.uk/ena/data/view/%s&set=true&display=fasta" % prefix

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

# def paired_fastq_url(acc):
#     """
#     Generate a URL for an SRA accession at EBI.
#     """
#     base = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq'
#     subdir = "/00%s" % acc[-1:] if len(acc) == 10 else ''
#     url = "%s/%s%s/%s/%s" % ( base, acc[:6], subdir, acc, acc )
#     return ["%s_%d.fastq.gz" % (url,i) for i in range(1,3)]
#
# def fastq_url(acc):
#     """
#     Generate a URL for an SRA accession at EBI.
#     """
#     # TODO: amend rule to allow length > 10
#     base = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq'
#     subdir = "/00%s" % acc[-1:] if len(acc) == 10 else ''
#     return "%s/%s%s/%s/%s.fastq.gz" % ( base, acc[:6], subdir, acc, acc )

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
        sra='SRR\d+'
    conda:
         '../envs/py3.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'enaDataGet -f fastq {params.sra}'
