rule fetch_assembly:
    """
    Fetch a remote assembly from EBI or NCBI.
    """
    output:
        fa='{assembly}.fasta'
    params:
        url=lambda wc: prepare_ncbi_assembly_url(config['assembly']['accession'],config['assembly']['alias'])
    wildcard_constraints:
        assembly='\w+'
    conda:
         '../envs/fetch.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'curl -s {params.url} | \
        pigz -d > {output.fa}'

# rule fetch_sra:
#     """
#     Fetch sra archive from EBI.
#     """
#     output:
#         temp('{sra}')
#     params:
#         sra = lambda wc: prepare_ebi_sra_url(wc.sra)
#     wildcard_constraints:
#         sra='\wRR\d+'
#     conda:
#          '../envs/fetch.yaml'
#     threads: 1
#     resources:
#         download=1,
#         threads=1
#     shell:
#         'wget -q {params.sra}'


rule fetch_fastq:
    """
    Fetch fastq file from EBI using aria2.
    """
    output:
        temp('{sra}{suff}.gz')
    params:
        url = lambda wc: prepare_ebi_sra_url(wc.sra,"%s.gz" % wc.suff)
    wildcard_constraints:
        sra='\wRR\d+',
        suff='[_\dsubread]*\.fastq'
    conda:
         '../envs/fetch.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        'aria2c -c --max-connection-per-server=8 --min-split-size=1M -s 8 {params.url}'

rule subsample_fastq:
    """
    Subsample large fastq files to reduce mapping time.
    """
    input:
        fastq=lambda wc: list_read_files(wc.sra,reads,False),
    output:
        temp("{sra}{suff}.gz")
    params:
        cmd = lambda wc: generate_subsample_command(wc.sra,reads)
    wildcard_constraints:
        sra='\wRR\d+',
        suff='[_\dsubread]*\.subsampled.fastq'
    conda:
         '../envs/blast.yaml'
    threads: 1
    resources:
        download=1,
        threads=1
    shell:
        '{params.cmd[0]} {input.fastq[0]} {params.cmd[1]} {output}'


# rule extract_paired_fastq:
#     """
#     Extract paired fastq files from sra archive.
#     """
#     input:
#         temp('{sra}')
#     output:
#         temp('{sra}_1.fastq.gz'),
#         temp('{sra}_2.fastq.gz')
#     conda:
#          '../envs/fetch.yaml'
#     threads: 1
#     resources:
#         download=1,
#         threads=1
#     shell:
#         'fastq-dump --split-files {input} && \
#         pigz {input}_*.fastq'
#
# rule extract_fastq:
#     """
#     extract fastq files from sra archive.
#     """
#     input:
#         temp('{sra}')
#     output:
#         temp('{sra}.fastq.gz')
#     wildcard_constraints:
#         sra='\wRR\d+'
#     conda:
#          '../envs/fetch.yaml'
#     threads: 1
#     resources:
#         download=1,
#         threads=1
#     shell:
#         'fastq-dump {input} && \
#         pigz {input}.fastq'
