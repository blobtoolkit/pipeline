rule unchunk_blast:
    """
    Unchunk chunked blast results.
    """
    input:
        '{assembly}.{type}.{name}.root.{root}{masked}.out.raw'
    output:
        '{assembly}.{type}.{name}.root.{root}{masked}.out'
    wildcard_constraints:
        root = r'\d+',
        type = '(blastn|diamond)',
        masked = r'.[fm][ulins\d\.]+',
        assembly = r'\w+'
    params:
        max_target_seqs = lambda wc: similarity[wc.name]['max_target_seqs']
    conda:
        '../envs/py3.yaml'
    threads: get_threads('unchunk_blast', 1)
    log:
        'logs/{assembly}/unchunk_blast/{type}.{name}.root.{root}{masked}.log'
    benchmark:
        'logs/{assembly}/unchunk_blast/{type}.{name}.root.{root}{masked}.benchmark.txt'
    script:
        '../scripts/unchunk_blast.py'
