#!/usr/bin/env python3

import math
import sys
import time
import shlex
import logging
from random import shuffle
from itertools import groupby
from multiprocessing import Pool
from subprocess import Popen, PIPE, run

logging.basicConfig(filename=snakemake.log[0], level=logging.WARNING)

try:
    FASTAFILE = snakemake.input.fasta
    BLAST_DB = snakemake.params.db
    TAXIDS = snakemake.input.taxids
    CHUNK = int(snakemake.params.chunk)
    OVERLAP = int(snakemake.params.overlap)
    MAXCHUNKS = int(snakemake.params.max_chunks)
    THREADS = int(snakemake.threads)
    EVALUE = str(snakemake.params.evalue)
    TARGETS = int(snakemake.params.max_target_seqs)
    PATH = snakemake.params.path
    OUTFILE = snakemake.output.raw
    NOHIT = snakemake.output.nohit
except Exception as err:
    logger.error(err)
    exit(1)

def chunk_size(value):
    """Calculate nice value for chunk size."""
    mag = math.floor(math.log10(value))
    first = int(str(value)[:2]) + 1
    chunk_size = first * pow(10, mag-2)
    return chunk_size


def split_list(input, size):
    """Yield successive subsets from list."""
    for i in range(0, len(input), size):
        yield input[i:i + size]


def chunk_fasta(fastafile, chunk=math.inf, overlap=0, max_chunks=math.inf):
    """
    Read FASTA file one sequence at a time and split long sequences into chunks
    """
    cmd = "cat %s" % fastafile
    # TODO: read gzipped files if needed
    # cmd = "pigz -dc %s" % fastafile
    title = ''
    seq = ''
    segment = chunk + overlap
    with Popen(shlex.split(cmd), encoding='utf-8', stdout=PIPE, bufsize=4096) as proc:
        faiter = (x[1] for x in groupby(proc.stdout, lambda line: line[0] == '>'))
        for header in faiter:
            title = header.__next__()[1:].strip().split()[0]
            seq = ''.join(map(lambda s: s.strip(),faiter.__next__()))
            seq_length = len(seq)
            if seq_length > segment:
                n = (seq_length + chunk) // chunk
                if n > max_chunks:
                    chunk = chunk_size(seq_length, max_chunks)
                    n = max_chunks
                for i in range(0, seq_length, chunk):
                    subseq = seq[i:i+segment]
                    yield {'title':title,'seq':subseq,'chunks':n,'start':i}
            else:
                yield {'title':title,'seq':seq,'chunks':1,'start':0}


def run_blast(seqs,db,evalue='1e-25',targets=10):
    """Run blast on seqs."""
    cmd = PATH+"/blastn -db %s \
              -outfmt \"6 qseqid staxids bitscore std\" \
              -max_target_seqs %d \
              -max_hsps 1 \
              -evalue %s \
              -taxidlist %s" % (db,targets,evalue,TAXIDS)
    input = ''
    for seq in seqs:
        input += ">%s_-_%d\n" % (seq['title'],seq['start'])
        input += "%s\n" % seq['seq']
    p = run(shlex.split(cmd), stdout=PIPE, stderr=PIPE, input=input, encoding='ascii')
    return p


if __name__ == '__main__':
    try:
        seqs = []
        names = set()
        for seq in chunk_fasta(FASTAFILE,CHUNK,OVERLAP):
            # ofh.write(">%s_-_%d\n" % (seq['title'],seq['start']))
            # ofh.write("%s\n" % seq['seq'])
            names.add(seq['title'])
            seqs.append((seq))
        n_chunks = len(seqs)
        shuffle(seqs)
        subset_length = math.ceil(n_chunks / THREADS)
        min_length = subset_length / THREADS
        while subset_length > THREADS and subset_length > min_length:
           subset_length //= 2
        pool = Pool(THREADS)
        jobs = []
        output = []

        def blast_callback(p):
            global output
            result = ''
            for line in p.stdout.strip('\n').split('\n'):
                fields = line.split('\t')
                if fields[0]:
                    output.append('\t'.join(fields))

        for subset in split_list(seqs, subset_length):
            proc = pool.apply_async(run_blast, (subset,BLAST_DB), callback=blast_callback)
            jobs.append(proc)
        pool.close()
        pool.join()
        for job in jobs:
            job.wait()
        with open(OUTFILE, 'w') as ofh:
            ofh.writelines('\n'.join(output))
        for line in output:
            name = line.split('_-_')[0]
            if name in names:
                names.remove(name)
        with open(NOHIT, 'w') as ofh:
            ofh.writelines('\n'.join(names))
    except Exception as err:
        logger.error(err)
        exit(1)
