#!/usr/bin/env python3

import subprocess
import os
import shutil
import gzip
import sys
import re
import shlex
from subprocess import Popen, run, PIPE
from itertools import groupby
from collections import defaultdict
from multiprocessing import Pool

# # set variables from snakemake params, wildcards, input and threads
# CHUNK = snakemake.params.chunk
# TMPDIR = "%s/%s/" % (snakemake.params.tmpdir,snakemake.wildcards.name)
# PATH = snakemake.wildcards.path
# DBTITLE = snakemake.wildcards.name
# THREADS = snakemake.threads
# FASTAFILE = snakemake.input.fa
# MAPFILE = snakemake.input.idmap
# OUTDIR = "%s/split/%s" % (PATH,NAME)

# set variables from snakemake params, wildcards, input and threads
PATH = '/ceph/software/databases/ncbi_2018_02/split/test'
DBTITLE = 'test.root.1.masked.28641'
THREADS = 4
TOOL = 'blast'
DBTYPE = 'nucl'
ROOT = str(1)
MASKS = [28641]
NODES = '/ceph/software/databases/ncbi_taxonomy/nodes.dmp'
TMPDIR = 'test'
OUTDIR = "%s/" % (DBTITLE)

if os.path.exists(TMPDIR):
    createtmp = False
else:
    createtmp = True
    os.makedirs(os.path.dirname("%s/" % TMPDIR), exist_ok=True)

TAXIDS = set()

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
    return 1

def mask_accessions(subset):
    """
    filter lines from idmap that are not in set of masked taxids,
    return list of accessions
    """
    accessions = []
    mapfile = "%s/%s.taxid_map.gz" % (PATH,subset)
    with Popen(['pigz', '-dc', mapfile], stdout=PIPE, encoding='utf-8', bufsize=4096) as proc:
        for line in proc.stdout:
            l = line.rstrip('\n').split()
            if len(l) > 1 and l[1] in TAXIDS:
                accessions.append(l[0])
    return accessions

def get_filename(file_):
    return "/dev/fd/{}".format(file_.fileno())

def get_stdout_fds(*processes):
    return tuple(p.stdout.fileno() for p in processes)

def make_makeblastbd_cmd(seqfile,dbtype,outdir,title,mapfile):
    cmd = ("makeblastdb -in %s -dbtype %s -title %s -out %s/%s -parse_seqids -taxid_map %s" %
            (seqfile,dbtype,title,outdir,title,mapfile))
    return shlex.split(cmd)

def makeblastdb(subset):
    print(subset)
    listfile = "%s/%s_%s.list" % (TMPDIR,DBTITLE,subset)
    with open(listfile,'w') as fh:
        fh.write('\n'.join(mask_accessions(subset))+'\n')
    seqfile = "%s/%s.fa.gz" % (PATH,subset)
    # seqtkpipe = Popen(
    #     ['seqtk','subseq',seqfile,listfile],
    #     stdout=PIPE,encoding='utf-8')
    # seqtkfile = get_filename(seqtkpipe.stdout)
    with open("%s.fa" % listfile,'w') as fh:
        seqtkpipe = Popen(
            ['seqtk','subseq',seqfile,listfile],
            stdout=fh,encoding='utf-8')
        print(seqtkpipe.stdout)
    # seqtkfile = get_filename(seqtkpipe.stdout)

    return
    mapfile = "%s/%s.taxid_map.gz" % (PATH,subset)
    idpipe = Popen(
        ['pigz','-dc',mapfile],
        stdout=PIPE,encoding='utf-8')
    idfile = get_filename(idpipe.stdout)
    cmd = make_makeblastbd_cmd(seqtkfile,DBTYPE,OUTDIR,"%s_%s" % (DBTITLE,subset),idfile)
    blast = run(
        cmd,stdout=PIPE,encoding='utf-8',
        pass_fds=get_stdout_fds(idpipe,seqtkpipe))
    os.remove(listfile)
    return blast.returncode

def main():
    graph = node_graph(NODES)
    masked = graph_to_masked_taxids(graph,ROOT,MASKS)
    os.makedirs(os.path.dirname(OUTDIR), exist_ok=True)
    subsets = [str(i).zfill(2) for i in range(0,2)]
    with Pool(THREADS) as p:
        p.map(makeblastdb,subsets)

if __name__ == '__main__':
    main()
