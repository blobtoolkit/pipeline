rule make_taxid_list:
    """
    Generate a list of taxids containing all descendants of a specified root,
    optionally with one or more lineages masked.
    """
    input:
        nodes = "%s/nodes.dmp" % config['settings']['taxonomy']
    output:
        '{name}.root.{root}{masked}.taxids',
        '{name}.root.{root}{masked}.negative.taxids'
    wildcard_constraints:
        root = r'\d+'
    params:
        mask_ids = lambda wc: similarity[wc.name]['mask_ids'],
        db = lambda wc: str("%s.root.%s%s" % (wc.name, wc.root, wc.masked)),
        taxdump = taxdump_dir
    conda:
        '../envs/py3.yaml'
    threads: get_threads('make_taxid_list', 1)
    log:
        'logs/make_taxid_list/{name}.root.{root}{masked}.log'
    benchmark:
        'logs/make_taxid_list/{name}.root.{root}{masked}.benchmark.txt'
    script:
        '../scripts/make_taxid_list.py'
