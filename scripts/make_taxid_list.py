#!/usr/bin/env python3

import os
import re
from collections import defaultdict

# set variables from snakemake params, wildcards, input and threads
DBTITLE = snakemake.params.db
ROOT = int(snakemake.wildcards.root)
MASKS = set(snakemake.params.mask_ids + [32630, 111789, 6])  # mask synthetic constructs by default
NODES = snakemake.input.nodes


def node_graph(nodes_file):
    """
    Read an NCBI nodes.dmp file and return a dict of child taxids and ranks for
    each taxid with descendants.
    """
    graph = defaultdict(dict)
    if os.path.isfile(nodes_file):
        with open(nodes_file, 'r') as fh:
            lines = fh.readlines()
            for l in lines:
                parts = re.split(r'\t\|\t', l)
                tax_id = int(parts[0])
                parent_id = int(parts[1])
                rank = parts[2]
                if parent_id == tax_id:
                    continue
                graph[parent_id].update({tax_id: rank})
    return graph


def graph_to_masked_taxids(graph, root, masks):
    """
    Generate a set containing masked taxids from a root.
    """
    taxids = set()

    def descend(root):
        """
        Iteratively descend from a root to generate a set of taxids
        unless the child taxid is in the list of taxids to mask.
        """
        if root in graph:
            for child, rank in graph[root].items():
                if masks and int(child) in masks:
                    continue
                taxids.add(child)
                descend(child)
    descend(root)

    return taxids


def main():
    masks = MASKS
    graph = node_graph(NODES)
    taxids = graph_to_masked_taxids(graph, ROOT, masks)
    listfile = "%s.taxids" % DBTITLE
    with open(listfile, 'w') as fh:
        fh.write("\n".join(str(taxid) for taxid in taxids))

    negative_ids = MASKS
    if ROOT != 1:
        negative_ids = negative_ids | graph_to_masked_taxids(graph, 1, set([ROOT]))
    for mask_id in masks:
        negative_ids = negative_ids | graph_to_masked_taxids(graph, mask_id, set())
    negative_listfile = "%s.negative.taxids" % DBTITLE
    with open(negative_listfile, 'w') as fh:
        fh.write("\n".join(str(taxid) for taxid in negative_ids))


if __name__ == '__main__':
    main()
