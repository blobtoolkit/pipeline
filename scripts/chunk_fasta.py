#!/usr/bin/env python3

import math
from itertools import groupby
from multiprocessing import Pool
from shlex import split
from subprocess import Popen, PIPE, run

FASTAFILE = snakemake.input
CHUNK = int(snakemake.params.chunk)
OVERLAP = int(snakemake.params.overlap)
OUTFILE = snakemake.output

def chunk_fasta(fastafile,chunk=math.inf,overlap=0):
    """
    Read FASTA file one sequence at a time and split long sequences into chunks
    """
    cmd = "cat %s" % fastafile
    # TODO: read gzipped files if needed
    # cmd = "pigz -dc %s" % fastafile
    title = ''
    seq = ''
    segment = chunk + overlap
    with Popen(split(cmd), encoding='utf-8', stdout=PIPE, bufsize=4096) as proc:
        faiter = (x[1] for x in groupby(proc.stdout, lambda line: line[0] == '>'))
        for header in faiter:
            title = header.__next__()[1:].strip().split()[0]
            seq = ''.join(map(lambda s: s.strip(),faiter.__next__()))
            seq_length = len(seq)
            if seq_length > segment:
                n = (seq_length + chunk) // chunk
                for i in range(0, seq_length, chunk):
                    subseq = seq[i:i+segment]
                    yield {'title':title,'seq':subseq,'chunks':n,'start':i}
            else:
                yield {'title':title,'seq':seq,'chunks':1,'start':0}

if __name__ == '__main__':
    with open(OUTPUT, 'w') as ofh:
        for seq in chunk_fasta(FASTAFILE,CHUNK,OVERLAP):
            ofh.write(">%s_-_%d\n" % (seq['title'],seq['start']))
            ofh.write("%s\n" % seq['seq'])
