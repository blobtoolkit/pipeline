#!/usr/bin/env python3
"""Update coordinates and remove seq name suffix from chunked blast results."""

import re
import logging
from collections import defaultdict
from docopt import docopt

docs = """
Unchunk BLAST results.

Usage: ./unchunk_blast.py [--count INT] [--in TSV] [--out TSV]

Options:
    --in TSV     input filename.
    --out TSV    output filename.
    --count INT  number of results to keep per chunk. [Default: 10]
"""

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

if __name__ == '__main__':
    args = docopt(docs)
    try:
        args['--in'] = snakemake.input[0]
        args['--count'] = snakemake.params.max_target_seqs
        args['--out'] = snakemake.output[0]
    except NameError as err:
        logger.info(err)
        logger.info('Parsing parameters from command line')

    try:
        lines = defaultdict(dict)
        chunk_counts = defaultdict(int)
        with open(args['--in'], 'r') as fh:
            for line in fh.readlines():
                fields = line.split('\t')
                if fields[0]:
                    name, start = re.split('_-_', fields[0])
                    fields[0] = name
                    fields[3] = name
                    fields[9] = str(int(fields[9])+int(start))
                    fields[10] = str(int(fields[10])+int(start))
                    if start not in lines[name]:
                        lines[name][start] = []
                        chunk_counts[name] += 1
                    lines[name][start].append('\t'.join(fields))
        with open(args['--out'], 'w') as ofh:
            for name in chunk_counts.keys():
                length = len(lines[name])
                n = int(args['--count'])
                for start in lines[name].keys():
                    for i in range(n):
                        if i < len(lines[name][start]):
                            ofh.write("%s" % (lines[name][start][i]))
    except Exception as err:
        logger.error(err)
        exit(1)
