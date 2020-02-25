rule make_masked_lists:
    """
    Generate a list of accessions needed to create a custom
    database containing all descendants of a specified root, optionally
    with one or more lineages masked.
    """
    input:
        split = lambda wc: "%s/split/%s.done" % (similarity[wc.name]['local'], wc.name),
        taxids= '{name}.root.{root}{masked}.negative.taxids'
    output:
        'blast/{name}.root.{root}{masked}.lists'
    wildcard_constraints:
        root = r'\d+'
    params:
        db = lambda wc: str("%s.root.%s%s" % (wc.name, wc.root, wc.masked)),
        indir = lambda wc: "%s/split/%s" % (similarity[wc.name]['local'], wc.name),
        chunk = config['settings']['chunk']
    conda:
        '../envs/py3.yaml'
    threads: get_threads('make_masked_lists', maxcore)
    log:
        'logs/make_masked_lists/{name}.root.{root}{masked}.log'
    benchmark:
        'logs/make_masked_lists/{name}.root.{root}{masked}.benchmark.txt'
    resources:
        threads = get_threads('make_masked_lists', maxcore)
    script:
        '../scripts/make_masked_lists.py'
