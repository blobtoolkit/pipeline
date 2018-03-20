#!/usr/bin/env python3

import os
import gzip
import pandas as pd
from collections import defaultdict
from pathlib import Path

# set variables from snakemake params, wildcards, input and threads
ROOT = int(snakemake.wildcards.root)
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
        if root in graph:
            for child,rank in graph[root].items():
                if masks and child in masks:
                    continue
                taxid = str(child)
                taxid_file = "%s/%s/%s" %(indir,taxid[-2:].zfill(2),taxid)
                if Path("%s.fa" % taxid_file).is_file():
                    filenames[taxid[-2:]].append(taxid_file)
                descend(child)
        return

    descend(root)
    return filenames

def node_graph(nodes_file):
    """
    Read an NCBI nodes.dmp file and return a dict of child taxids and ranks for
    each taxid with descendants.
    """
    df = pd.read_table(nodes_file,sep='\t\|\t',header=None,engine='python')
    df.rename(columns={0:'tax_id',1:'parent_id',2:'rank'},inplace=True)
    graph = defaultdict(dict)
    for row in df.itertuples():
        tax_id = getattr(row,'tax_id')
        parent_id = getattr(row,'parent_id')
        rank = getattr(row,'rank')
        if parent_id == tax_id:
            continue
        graph[parent_id].update({tax_id:rank})
    return graph

# Write a file containing a list of per-taxon sequence filenames needed to
# create a custom database containing all descendants of a specified root,
# optionally with one or more lineages masked.
graph = node_graph(NODES)
filenames = graph_to_filenames(graph,ROOT,MASKS,INDIR,PREFIX)
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
