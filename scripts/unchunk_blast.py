#!/usr/bin/env python3

import re
from collections import defaultdict

INFILE = snakemake.input
OUTFILE = snakemake.output
COUNT = snakemake.params.max_target_seqs

if __name__ == '__main__':
    lines = defaultdict(dict)
    chunk_counts = defaultdict(int)
    with open(INFILE, 'r') as fh:
        for line in fh.readlines():
            fields = line.split('\t')
            if fields[0]:
                name,start = re.split('_-_', fields[0])
                fields[0] = name
                fields[3] = name
                fields[9] = str(int(fields[9])+int(start))
                fields[10] = str(int(fields[10])+int(start))
                if start not in lines[name]:
                    lines[name][start] = []
                    chunk_counts[name] += 1
                lines[name][start].append('\t'.join(fields))
    with open(OUTFILE, 'w') as ofh:
        for name in chunk_counts.keys():
            l = len(lines[name])
            n = (COUNT+l-1)// l
            for start in lines[name].keys():
                for i in range(n):
                    if i < len(lines[name][start]):
                        ofh.write("%s" % (lines[name][start][i]))
