#!/usr/bin/env python3

import os
import re

INFILES = snakemake.input
OUTFILE = str(snakemake.output)

total = mapped = unmapped = 0
blobtools = source = ''
coverage = {}
order = []

for file in INFILES:
    with open(file,'r') as fh:
        for line in fh:
            line = line.lstrip('\n')
            if line.startswith('#'):
                if line.startswith('## blobtools'):
                    blobtools = re.search(r'(v.+)',line).group(1)
                elif line.startswith('## Total Reads'):
                    total += int(re.search(r'(\d+)',line).group(1))
                elif line.startswith('## Mapped Reads'):
                    mapped += int(re.search(r'(\d+)',line).group(1))
                elif line.startswith('## Unmapped Reads'):
                    unmapped += int(re.search(r'(\d+)',line).group(1))
                elif line.startswith('## Source'):
                    source += re.search(r':(.+)',line).group(1)
            else:
                contig_id,read_cov,base_cov = line.split()
                if contig_id not in coverage:
                    coverage[contig_id] = {'read_cov':int(read_cov),'base_cov':float(base_cov)}
                    order.append(contig_id)
                else:
                    coverage[contig_id]['read_cov'] += int(read_cov)
                    coverage[contig_id]['base_cov'] += float(base_cov)

header = """## blobtools %s
## Total Reads = %d
## Mapped Reads = %d
## Unmapped Reads = %d
## Source(s) : %s
# contig_id	read_cov	base_cov
""" % (blobtools,total,mapped,unmapped,source)

with open(OUTFILE,'w') as fh:
    fh.write(header)
    for contig_id in order:
        fh.write("%s\t%d\t%.3f\n" % (contig_id,coverage[contig_id]['read_cov'],coverage[contig_id]['base_cov']))
