rule run_blastn:
    """
    Run NCBI blastn to search nucleotide database with assembly query.
    """
    input:
        # fasta = '{assembly}.windowmasker.fasta',
        fasta = '{assembly}.chunks.fasta',
        db = lambda wc: "%s/%s.nal" % (similarity[wc.name]['local'], wc.name),
        taxids = '{name}.root.{root}{masked}.taxids'
    output:
        nohit = '{assembly}.blastn.{name}.root.{root}{masked}.nohit' if keep else temp('{assembly}.blastn.{name}.root.{root}{masked}.nohit'),
        raw = '{assembly}.blastn.{name}.root.{root}{masked}.out.raw'
    wildcard_constraints:
        root = '\d+',
        masked = '.[fm][ulins\d\.]+'
    params:
        dir = ncbi_dir,
        evalue = lambda wc: similarity[wc.name]['evalue'],
        max_target_seqs = lambda wc: similarity[wc.name]['max_target_seqs'],
        multiprocessing = True if 'multiprocessing' in config['settings'] else False,
        chunk = 0,
        overlap = config['settings']['blast_overlap'],
        max_chunks = config['settings']['blast_max_chunks']
    conda:
        '../envs/pyblast.yaml'
    threads: get_threads('run_blastn', maxcore)
    log:
        'logs/{assembly}/run_blastn/{name}.root.{root}{masked}.log'
    benchmark:
        'logs/{assembly}/run_blastn/{name}.root.{root}{masked}.benchmark.txt'
    script:
        '../scripts/blast_wrapper.py'
