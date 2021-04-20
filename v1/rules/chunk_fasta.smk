rule chunk_fasta:
    """
    Split long contigs into chunks.
    """
    input:
        fasta = '{assembly}.windowmasker.fasta'
    output:
        '{assembly}.chunks.fasta' if keep else temp('{assembly}.chunks.fasta')
    params:
        chunk = config['settings']['blast_chunk'],
        overlap = config['settings']['blast_overlap'],
        max_chunks = config['settings']['blast_max_chunks']
    conda:
        '../envs/blobtools2.yaml'
    threads: get_threads('chunk_fasta', 1)
    log:
        'logs/{assembly}/chunk_fasta.log'
    benchmark:
        'logs/{assembly}/chunk_fasta.benchmark.txt'
    script:
        '../scripts/chunk_fasta.py'
