#!/usr/bin/env python3

import json
import math
import time
from itertools import groupby
from multiprocessing import Pool
from shlex import split
from subprocess import Popen, PIPE, run

FASTAFILE = snakemake.input.fasta
BLAST_DB = snakemake.params.db
CHUNK = int(snakemake.params.blast_chunk)
OVERLAP = int(snakemake.params.blast_overlap)
THREADS = int(snakemake.threads)
TARGETS = int(snakemake.params.max_target_seqs)
EVALUE = str(snakemake.params.evalue)
TAXIDS = snakemake.input.taxids
PATH = snakemake.params.path
OUTFILE = snakemake.output[0]

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

def run_blast(seq,db='nt',evalue='1e-25',targets=10):
    """
    run blast on single (chunked) sequences
    """
    cmd = "%s/blastn -db %s \
              -outfmt \"6 qseqid staxids bitscore std\" \
              -max_target_seqs %d \
              -max_hsps 1 \
              -evalue %s \
              -taxidlist %s" % (PATH,db,targets,evalue,TAXIDS)
    p = run(split(cmd), stdout=PIPE, stderr=PIPE, input=">%s\n%s\n" % (seq['title'],seq['seq']), encoding='ascii')
    output = []
    result = ''
    print(p.stderr)
    for line in p.stdout.strip('\n').split('\n'):
        fields = line.split('\t')
        if fields[0]:
            fields[9] = str(int(fields[9])+seq['start'])
            fields[10] = str(int(fields[10])+seq['start'])
            output.append('\t'.join(fields))
            result = "%s\n" % '\n'.join(output)
    return result

if __name__ == '__main__':
    p = Pool(THREADS)
    results = []
    for seq in chunk_fasta(FASTAFILE,CHUNK,OVERLAP):
        targets = max(TARGETS // seq['chunks'],1)
        result = p.apply_async(run_blast,(seq,),{'db':BLAST_DB,'evalue':EVALUE,'targets':targets})
        results.append(result)
    p.close()
    p.join()
    results = ''.join([result.get() for result in results])
    with open(OUTFILE, 'w') as ofh:
        ofh.write(results)
