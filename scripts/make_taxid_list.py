#!/usr/bin/env python3

import os
import re
from collections import defaultdict

# set variables from snakemake params, wildcards, input and threads
DBTITLE = snakemake.params.db
ROOT = str(snakemake.wildcards.root)
MASKS = snakemake.params.mask_ids + [32630,111789,6] # mask synthetic constructs by default
NODES = snakemake.input.nodes

# mask synthetic constructs
TAXIDS = {'32630','111789','6'}

def node_graph(nodes_file):
    """
    Read an NCBI nodes.dmp file and return a dict of child taxids and ranks for
    each taxid with descendants.
    """
    graph = defaultdict(dict)
    if os.path.isfile(nodes_file):
        with open(nodes_file,'r') as fh:
            lines = fh.readlines()
            for l in lines:
                parts = re.split(r'\t\|\t',l)
                tax_id = parts[0]
                parent_id = parts[1]
                rank = parts[2]
                if parent_id == tax_id:
                    continue
                graph[parent_id].update({tax_id:rank})
    return graph

def graph_to_masked_taxids(graph,root,masks):
    """
    Generate a set containing masked taxids from a root.
    """
    global TAXIDS

    def descend(root):
        """
        Iteratively descend from a root to generate a set of taxids
        unless the child taxid is in the list of taxids to mask.
        """
        if root in graph:
            for child,rank in graph[root].items():
                if masks and int(child) in masks:
                    continue
                TAXIDS.add(child)
                descend(child)
    descend(root)

def graph_to_masked_taxids(graph,root,masks):
    """
    Generate a set containing masked taxids from a root.
    """
    global TAXIDS

    def descend(root):
        """
        Iteratively descend from a root to generate a set of taxids
        unless the child taxid is in the list of taxids to mask.
        """
        if root in graph:
            for child,rank in graph[root].items():
                if masks and int(child) in masks:
                    continue
                TAXIDS.add(child)
                descend(child)
    descend(root)

def main():
    graph_to_masked_taxids(node_graph(NODES),ROOT,MASKS)
    listfile = "%s.taxids" % DBTITLE
    with open(listfile, 'w') as fh:
        fh.write("\n".join(TAXIDS))

if __name__ == '__main__':
    main()
