#!/usr/bin/env python3

import os
import shutil
import gzip
import sys
import re
from itertools import groupby
from collections import defaultdict
from subprocess import Popen, call, PIPE
from multiprocessing import Pool
import timeit

# set variables from snakemake params, wildcards, input and threads
CHUNK = snakemake.params.chunk
NAME = snakemake.wildcards.name
TMPDIR = "%s/%s" % (snakemake.params.tmpdir,NAME)
THREADS = snakemake.threads
FASTAFILE = snakemake.input.fa
MAPFILE = snakemake.input.idmap
OUTDIR = snakemake.params.outdir

# CHUNK = 1000000
# PATH = '/ceph/software/databases/ncbi_2018_02'
# NAME = 'nr'
# THREADS = 8
# FASTAFILE = '/ceph/software/databases/ncbi_2018_02/full/nr.fa.gz'
# MAPFILE = '/ceph/software/databases/ncbi_2018_02/full/nr.taxid_map.gz'
# TMPDIR = '/dev/shm/nr'
# OUTDIR = '/scratch/rchallis/nr'

# create temporary directory to write files into
if os.path.exists(TMPDIR):
    createtmp = False
else:
    createtmp = True
    os.makedirs(os.path.dirname("%s/" % TMPDIR), exist_ok=True)

# read the id mapping file into a dict
mapping = {}
with Popen(['pigz', '-dc', MAPFILE], stdout=PIPE, encoding='utf-8', bufsize=4096) as proc:
    for line in proc.stdout:
        l = line[:-1].split()
        mapping.update({l[0]:l[1]})

def fasta_chunks(fastafile,chunk):
    """
    Read batches of sequences from a gzipped fasta file
    return dict of id mapping and fasta sequence by taxid
    """
    by_taxid_suffix = defaultdict(lambda: {'acc':[],'fa':[]})
    ctr = 0
    with Popen(['pigz', '-dc', fastafile], stdout=PIPE, encoding='utf-8', bufsize=4096) as proc:
        faiter = (x[1] for x in groupby(proc.stdout, lambda line: line[0] == '>'))
        for header in faiter:
            taxids = header_to_taxid(header.__next__()[1:].strip())
            seq = ''.join(map(lambda s: s.strip(),faiter.__next__()))
            if taxids:
                ctr += 1
                for taxid,acc in taxids.items():
                    suffix = taxid[-2:].zfill(2)
                    by_taxid_suffix[suffix]['acc'].append("%s\t%s\n" % (acc,taxid))
                    by_taxid_suffix[suffix]['fa'].append(">%s\n%s\n" % (acc,seq))
            if ctr % CHUNK == 0:
                yield by_taxid_suffix
                by_taxid_suffix = defaultdict(lambda: {'acc':[],'fa':[]})
    return by_taxid_suffix

def header_to_taxid(header):
    """
    Split each sequence header on '\\x01' and return a dict with one fasta
    format sequence per taxid
    """
    seen_taxids = {}
    taxids = {}
    for entry in header.split('\x01'):
        acc = entry.split(' ',1)[0]
        try:
            taxid = mapping[acc]
        except:
            continue
            # failure may be due to mismatched naming convention with pdb
            # accessions, removing the last character can help but prone to
            # error so just ignore anything that doesn't match
            # if acc[:-1] in mapping:
            #     taxid = mapping[acc[:-1]]
            # else:
            #     # TODO: log missing entries
            #     continue
        if taxid not in seen_taxids:
            seen_taxids[taxid] = True
            taxids[taxid] = acc
    return taxids

by_tax_suffix = {}

def write_files(suffix):
    """
    Write a list of sequences to a gzipped fasta file and corresponding taxid
    map file
    """
    l = by_tax_suffix[suffix]
    fileroot = "%s/%s" % (TMPDIR,suffix)
    fafile = "%s.fa.gz" % fileroot
    with open(fafile, 'a') as fh:
        with Popen(['pigz', '-c'], stdin=PIPE, stdout=fh, encoding='utf-8', bufsize=4096) as proc:
            proc.stdin.write(''.join(l['fa']))
    accfile = "%s.taxid_map.gz" % fileroot
    with open(accfile, 'a') as fh:
        with Popen(['pigz', '-c'], stdin=PIPE, stdout=fh, encoding='utf-8', bufsize=4096) as proc:
            proc.stdin.write(''.join(l['acc']))
    return 1

def make_redundant(fastafile,chunk,threads,outdir,tmpdir):
    """
    read a non-redundant fasta file, duplicate sequences where necessary and
    write one sequence file per taxid suffix
    """
    global by_tax_suffix
    fiter = fasta_chunks(fastafile,chunk)
    ctr = 0
    start_time = timeit.default_timer()
    total_time = 0
    for by_taxid_suffix in fiter:
        ctr += 1
        elapsed = timeit.default_timer() - start_time
        print(" - prepared (%s)" % str(elapsed), file=sys.stderr)
        # write sequences and taxid maps to file
        by_tax_suffix = by_taxid_suffix
        with Pool(processes=threads) as p:
            p.map(write_files,by_taxid_suffix.keys())
        elapsed = timeit.default_timer() - start_time
        print(" - written (%s)" % str(elapsed), file=sys.stderr)
        total_time += elapsed
        print("Processed %d sequences in %s seconds" % (chunk * ctr,str(total_time)), file=sys.stderr)
        start_time = timeit.default_timer()
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
