def prepare_ncbi_assembly_url(accession,name):
    """
    Generate a URL for an assembly at NCBI.
    """
    base = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all'
    acc = accession.replace('_','').split('.',1)[0]
    path = '/'.join(acc[i:i+3] for i in range(0, len(acc), 3))
    asm = "%s_%s" % ( accession, name )
    url = "%s/%s/%s/%s_genomic.fna.gz" % ( base, path, asm, asm )
    return url

def prepare_ebi_assembly_url(insdc):
    """
    Generate a URL for an assembly at EBI.
    """
    base = 'ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public'
    url = "%s/%s/%s.fasta.gz" % ( base, insdc[0:2].lower(), insdc[0:6] )
    return url

def prepare_assembly_url(config):
    """
    Generate an appropriate URL for an assembly accession depending on source.
    """
    if config['assembly']['source'] == 'ebi':
        return prepare_ebi_assembly_url(config['assembly']['insdc'])
    elif config['assembly']['source'] == 'ncbi':
        return prepare_ncbi_assembly_url(config['assembly']['accession'],config['assembly']['name'])

rule fetch_assembly:
    """
    Fetch a remote assembly from EBI or NCBI.
    """
    output:
        temp('{assembly}.fna')
    params:
        url=prepare_assembly_url(config)
    threads: 1
    resources:
        download=1
    shell:
        'curl {params.url} | pigz -dc > {output}'

def prepare_ncbi_sra_url(acc):
    """
    Generate a URL for an SRA accession at NCBI.
    """
    base = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra'
    url = "%s/%s/%s/%s/%s.sra" % ( base, acc[:3], acc[:6], acc, acc )
    return url

def prepare_ebi_sra_url(acc):
    """
    Generate a URL for an SRA accession at EBI.
    """
    base = 'ftp://ftp.sra.ebi.ac.uk/vol1'
    subdir = "/00%s" % acc[-1:] if len(acc) == 10 else ''
    url = "%s/%s/%s%s/%s" % ( base, acc[:3].lower(), acc[:6], subdir, acc )
    return url

def prepare_sra_url(acc):
    """
    Generate an appropriate URL for an SRA accession depending on source.
    """
    if config['reads']['source'] == 'ebi':
        return prepare_ebi_sra_url(acc)
    elif config['reads']['source'] == 'ncbi':
        return prepare_ncbi_sra_url(acc)

rule fetch_interleaved_sra:
    """
    Fetch an SRA file from EBI or NCBI for paired data.
    """
    params:
        url=lambda wildcards:
            dict([i, prepare_sra_url(i)] for i in config['reads']['paired'])[wildcards.input]
    output:
        temp('{input}.interleaved.sra')
    threads: 1
    resources:
        download=1
    shell:
        'curl {params.url} > {output}'

rule interleaved_fastq_dump:
    """
    Extract an interleaved FASTQ file from an SRA dump.
    """
    input:
        '{sample}.interleaved.sra'
    output:
        temp('{sample}.interleaved.fastq')
    params:
        sample=lambda wc: wc.sample
    threads: 1
    shell:
        '{ENV} \
        fastq-dump {input}'
