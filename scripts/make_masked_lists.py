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
TAXID_FILE = snakemake.input.taxids
THREADS = snakemake.threads
OUTDIR = 'blast'


def read_taxids(filename):
    """
    read masked list of taxids from file
    """
    with open(filename, 'r') as fh:
        taxids = fh.read().splitlines()
    return taxids


def mask_accessions(subset, taxids):
    """
    filter lines from idmap that are not in set of masked taxids,
    return list of accessions
    """
    accessions = []
    mapfile = "%s/%s.taxid_map.gz" % (INDIR, subset)
    if os.path.isfile(mapfile):
        with Popen(['pigz', '-dc', mapfile], stdout=PIPE, encoding='utf-8', bufsize=4096) as proc:
            for line in proc.stdout:
                row = line.rstrip('\n').split()
                if len(row) > 1 and row[1] not in taxids:
                    accessions.append(row[0])
    return accessions


def write_lists(subset, taxids):
    listfile = "%s/%s_%s.accessions" % (OUTDIR, DBTITLE, subset)
    accessions = mask_accessions(subset, taxids)
    if accessions:
        with open(listfile, 'w') as fh:
            fh.write('\n'.join(accessions)+'\n')
        return subset
    return


def main():
    taxids = read_taxids(TAXID_FILE)
    os.makedirs(os.path.dirname("%s/" % OUTDIR), exist_ok=True)
    subsets = [(str(i).zfill(2), taxids) for i in range(0, 100)]
    with Pool(THREADS) as p:
        sets = list(p.starmap(write_lists, subsets))
    listfile = "%s/%s.lists" % (OUTDIR, DBTITLE)
    with open(listfile, 'w') as fh:
        for subset in sets:
            if subset:
                fh.write("%s\n" % subset)


if __name__ == '__main__':
    main()
