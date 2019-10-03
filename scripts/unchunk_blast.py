#!/usr/bin/env python3

import re
import logging
from collections import defaultdict

logger_config = {
    'level': logging.INFO,
    'format': '%(asctime)s [%(levelname)s] line %(lineno)d %(message)s',
    'filemode': 'w'
}
try:
    logger_config.update({'filename': snakemake.log[0]})
except NameError as err:
    pass
logging.basicConfig(**logger_config)
logger = logging.getLogger()

try:
    INFILE = snakemake.input[0]
    OUTFILE = snakemake.output[0]
    COUNT = snakemake.params.max_target_seqs
except Exception as err:
    logger.error(err)
    exit(1)

if __name__ == '__main__':
    try:
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
    except Exception as err:
        logger.error(err)
        exit(1)
