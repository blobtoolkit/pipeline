#!/usr/bin/env python3

import os
import re
from subprocess import Popen, PIPE
from collections import defaultdict
from multiprocessing import Pool

# set variables from snakemake params, wildcards, input and threads
CHUNK = snakemake.params.chunk
INDIR = snakemake.params.indir
DBTITLE = snakemake.params.db
ROOT = str(snakemake.wildcards.root)
MASKS = snakemake.params.mask_ids
NODES = snakemake.input.nodes
THREADS = snakemake.threads
OUTDIR = 'blast'

# # set variables from snakemake params, wildcards, input and threads
# INDIR = '/scratch/rchallis/nt'
# DBTITLE = 'test.root.1.masked.28641'
# THREADS = 8
# ROOT = str(1)
# MASKS = [28641]
# NODES = '/ceph/software/databases/ncbi_taxonomy/nodes.dmp'
# OUTDIR = 'blast'

TAXIDS = set()

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
    Generate a set containing masked of taxids from a root.
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

def mask_accessions(subset):
    """
    filter lines from idmap that are not in set of masked taxids,
    return list of accessions
    """
    accessions = []
    mapfile = "%s/%s.taxid_map.gz" % (INDIR,subset)
    if os.path.isfile(mapfile):
        with Popen(['pigz', '-dc', mapfile], stdout=PIPE, encoding='utf-8', bufsize=4096) as proc:
            for line in proc.stdout:
                l = line.rstrip('\n').split()
                if len(l) > 1 and l[1] in TAXIDS:
                    accessions.append(l[0])
    return accessions

def write_lists(subset):
    listfile = "%s/%s_%s.accessions" % (OUTDIR,DBTITLE,subset)
    accessions = mask_accessions(subset)
    if accessions:
        with open(listfile,'w') as fh:
            fh.write('\n'.join(accessions)+'\n')
        return subset
    return

def main():
    graph_to_masked_taxids(node_graph(NODES),ROOT,MASKS)
    os.makedirs(os.path.dirname("%s/" % OUTDIR), exist_ok=True)
    subsets = [str(i).zfill(2) for i in range(0,100)]
    with Pool(THREADS) as p:
        sets = list(p.map(write_lists,subsets))
    listfile = "%s/%s.lists" % (OUTDIR,DBTITLE)
    with open(listfile, 'w') as fh:
        for subset in sets:
            if subset:
                fh.write("%s\n" % subset)

if __name__ == '__main__':
    main()
