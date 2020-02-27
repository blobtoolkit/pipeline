#!/usr/bin/env python3

import re
import gzip
import os
import tarfile
import glob
import shutil

# Set variables from snakemake params/wildcards
TMPDIR = snakemake.params.tmpdir
TARFILE = snakemake.input[0]
NAME = snakemake.wildcards.name
OUTDIR = TARFILE.replace('.tar.gz', '')

# create temporary directory to write files into
if os.path.exists(TMPDIR):
    createtmp = False
else:
    createtmp = True
    os.makedirs(os.path.dirname(TMPDIR), exist_ok=True)

def regex_files(regex,members):
    """
    Function to filter a list of filenames to those matching a regex.
    """
    for tarinfo in members:
        if re.match(regex,tarinfo.name):
            yield tarinfo

# extract files from the uniprot tar archive
tar = tarfile.open(TARFILE)
# extract all protein fasta files (ignoring *.DNA_fasta* and *.additional_fasta*)
fasta_regex = r'.+\d\.fasta\.gz'
tar.extractall(path=TMPDIR,members=regex_files(fasta_regex,tar))
# extract all corresponding idmapping files
id_regex = r'.+\d\.idmapping\.gz'
tar.extractall(path=TMPDIR,members=regex_files(id_regex,tar))
tar.close()

# simplify fasta file headers and write to a single outfile
fasta_files = glob.glob("%s/*/*.fasta.gz" % TMPDIR)
for f in fasta_files:
    outlines = []
    with gzip.open(f, 'rt') as fh:
        inlines = fh.readlines()
        for l in inlines:
            if l[0] == '>':
                l = ">%s\n" % l.split('|')[1]
            outlines.append(l)
    with gzip.open("%s.fa.gz" % OUTDIR, 'at') as fh:
        fh.writelines(outlines)
    os.remove(f)

# extract all NCBI_TaxID lines from the idmapping files,
# reformat as BLAST taxid_map
# and write to a single outfile
mapping_files = glob.glob("%s/*/*.idmapping.gz" % TMPDIR)
for f in mapping_files:
    outlines = []
    regex = r'.+\tNCBI_TaxID\t.+'
    with gzip.open(f, 'rt') as fh:
        inlines = fh.readlines()
        for l in inlines:
            if re.match(regex,l):
                acc,tmp,taxid = l.strip().split('\t')
                outlines.append("%s\t%s\n" % (acc,taxid))
    with gzip.open("%s.taxid_map.gz" % OUTDIR, 'at') as fh:
        fh.writelines(outlines)
    os.remove(f)

# remove the temporary directory and any remaining content
if createtmp:
    shutil.rmtree(TMPDIR)
