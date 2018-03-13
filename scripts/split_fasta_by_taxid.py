#!/usr/bin/env python

# ===================================================================== #
# This script takes too long to run, particularly on nt.                #
# Need to improve efficiency but size of dataset caused problems with a #
# simple multiprocessing approach, hence commented out.                 #
# ===================================================================== #

# ===================================================================== #
# NB: runnning this script on nt will create directories with >100000   #
#     files. Would probably be best to split further to avoid           #
#     challenging filesystems (would require corresponding changes to   #
#     masked_list_by_root.py).                                          #
# ===================================================================== #

import os
import shutil
import gzip
import sys
from itertools import groupby
from collections import defaultdict
from subprocess import call
# from pathlib import Path
# from multiprocessing import Pool

# set variables from snakemake params, wildcards, input and threads
CHUNK = snakemake.params.chunk
# TODO: fix the unique tmpdir name as no snakemake.jobid
TMPDIR = "%s/%s/" % (snakemake.params.tmpdir,snakemake.jobid)
PATH = snakemake.wildcards.path
NAME = snakemake.wildcards.name
THREADS = snakemake.threads
FASTAFILE = snakemake.input.fa
MAPFILE = snakemake.input.idmap
OUTDIR = "%s/split/%s" % (PATH,NAME)

# create temporary directory to write files into
if os.path.exists(TMPDIR):
    createtmp = False
else:
    createtmp = True
    os.makedirs(os.path.dirname(TMPDIR), exist_ok=True)

# read the id mapping file into a dict
mapping = {}
ctr = 0;
with gzip.open(MAPFILE, 'rt', encoding='utf-8') as gz:
    for line in gz:
        if line.strip():
            ctr += 1
            a,b =  line.strip().split()
            mapping.update({a:b})

def fasta_chunks(fastafile,chunk):
    """
    Read batches of sequences from a gzipped fasta file
    return list of header/sequence tuples
    """
    fh = gzip.open(fastafile, 'rt', encoding='utf-8')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == '>'))
    seqs = []
    for header in faiter:
        header_str = header.__next__()[1:]
        seq = ''.join(map(lambda s: s.strip(),faiter.__next__()))
        seqs.append((header_str,seq))
        if len(seqs) % CHUNK == 0:
            yield seqs
            seqs = []
    return seqs

def duplicate_seqs(sequence):
    """
    Split each sequence header on '\\x01' and return a dict with one fasta
    format sequence per taxid
    """
    header,seq = sequence
    taxids = {}
    for entry in header.split('\x01'):
        acc = entry.split(' ',1)[0]
        try:
            taxid = mapping[acc]
        except:
            # failure may be due to mismatched naming convention with pdb
            # accessions, removing the last character can help
            if acc[:-1] in mapping:
                taxid = mapping[acc[:-1]]
            else:
                # TODO: log missing entries
                continue
        if taxid not in taxids:
            taxids[taxid] = ">%s\n%s\n" % (acc,seq)
    return taxids

def write_file(data):
    """
    Write a list of sequences to a gzipped fasta file
    """
    taxid,lines = data
    filepath = "%s/%s/%s.fa.gz" % (TMPDIR,taxid[-1],taxid)
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with gzip.open(filepath, 'at', encoding='utf-8') as gz:
        gz.writelines(lines)
    return 1

def make_redundant(fastafile,chunk,threads,outdir,tmpdir):
    """
    read a non-redundant fasta file, duplicate sequences where necessary and
    write one sequence file per taxid
    """
    fiter = fasta_chunks(fastafile,chunk)
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    if not tmpdir:
        tmpdir = outdir
    ctr = 0
    for seqs in fiter:
        ctr += 1
        print("processing %d sequences..." % (chunk * ctr), file=sys.stderr)
        by_taxid = defaultdict(list)
        for seq in seqs:
            taxids = duplicate_seqs(seq)
            for taxid,fasta in taxids.items():
                by_taxid[taxid].append(fasta)
        # TODO: Hit multiprocessing bug due to size of dataset when processing nt
        # need to find a way to reinstate multiprocessing as otherwise
        # processing nt can take ~ 3 days
        # os.makedirs(os.path.dirname(outdir), exist_ok=True)
        # with Pool(threads) as p:
        #     p.map(write_file,by_taxid.items())
        list(map(write_file,by_taxid.items()))
        print('\t...done', file=sys.stderr)
    if not tmpdir == outdir:
        # move files to outdir
        os.makedirs(os.path.dirname(outdir), exist_ok=True)
        call(['rsync','--remove-source-files','-a',"%s/" % tmpdir,"%s/" % outdir])
        # remove the temporary directory and any remaining content
        if createtmp:
            shutil.rmtree(tmpdir)

def main():
    make_redundant(FASTAFILE,CHUNK,THREADS,OUTDIR,TMPDIR)

if __name__ == '__main__':
    main()
