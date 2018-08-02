#!/usr/bin/env python3

import os
import gzip
import re
from collections import defaultdict
from pathlib import Path

# set variables from snakemake params, wildcards, input and threads
ROOT = str(snakemake.wildcards.root)
MASKS = snakemake.params.mask_ids
NODES = snakemake.input.nodes
INDIR = snakemake.params.indir
PREFIX = snakemake.wildcards.name

def graph_to_filenames(graph,root,masks,indir,prefix):
    """
    Generate a list of filenames for the sequences of nested taxids from a root.
    """
    filenames = defaultdict(list)

    def descend(root):
        """
        Iteratively descend from a root to generate a list of
        filenames unless the child taxid is in the list of taxids to mask.
        """
        print(root)
        if root in graph:
            print(root)
            for child,rank in graph[root].items():
                print(child)
                if masks and int(child) in masks:
                    continue
                taxid_file = "%s/%s/%s" %(indir,child[-2:].zfill(2),child)
                if Path("%s.fa" % taxid_file).is_file():
                    filenames[child[-2:]].append(taxid_file)
                descend(child)
        return

    descend(root)
    return filenames

def node_graph(nodes_file):
    """
    Read an NCBI nodes.dmp file and return a dict of child taxids and ranks for
    each taxid with descendants.
    """
    graph = defaultdict(dict)
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

# Write a file containing a list of per-taxon sequence filenames needed to
# create a custom database containing all descendants of a specified root,
# optionally with one or more lineages masked.
graph = node_graph(NODES)
filenames = graph_to_filenames(graph,ROOT,MASKS,INDIR,PREFIX)
if len(filenames.items()) > 0:
    suffix = ''
    if MASKS:
        suffix = ".minus.%s" % '.'.join(str(mask) for mask in MASKS)
    outdir = "%s.root.%s%s/" %(PREFIX,ROOT,suffix)
    os.makedirs(os.path.dirname(outdir), exist_ok=True)
    dir_list = "%s.root.%s%s/list" % (PREFIX,ROOT,suffix)
    with open(dir_list, 'w') as dl:
        for directory,files in filenames.items():
            file_list = "%s%s.list" %(outdir,directory)
            with open(file_list, 'w') as fh:
                fh.write('\n'.join(files))
            dl.write("%s\n" % file_list)
