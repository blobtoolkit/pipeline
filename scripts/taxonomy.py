#!/usr/bin/env python3

import os
import re
from collections import defaultdict

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

def parents_at_rank(graph,root,parent_rank):
    """
    loop through graph from root taxon, assigning leaf nodes to parent nodes at
    a given rank.
    """
    parents = {}
    def descend(root,parent):
        """
        Iteratively descend from a root to generate a set of taxids
        unless the child taxid is in the list of taxids to mask.
        """
        if root in graph:
            for child,rank in graph[root].items():
                if rank == parent_rank:
                    parent = child
                elif parent:
                    parents[child] = parent
                descend(child,parent)
    descend(str(root),None)
    return parents
