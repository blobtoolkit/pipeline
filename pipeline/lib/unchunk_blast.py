#!/usr/bin/env python3
"""
Unchunk BLAST results.

Update coordinates and remove seq name suffix from chunked blast results.

Usage: unchunk-blast [--count INT] --in TSV --out TSV

Options:
    --in TSV     input filename.
    --out TSV    output filename.
    --count INT  number of results to keep per chunk. [Default: 10]
"""

import logging
import re
import sys
from collections import defaultdict

from docopt import DocoptExit
from docopt import docopt

logger_config = {
    "level": logging.INFO,
    "format": "%(asctime)s [%(levelname)s] line %(lineno)d %(message)s",
    "filemode": "w",
}
try:
    logger_config.update({"filename": snakemake.log[0]})
except NameError as err:
    pass
logging.basicConfig(**logger_config)
logger = logging.getLogger()


def parse_args():
    """Parse snakemake args if available."""
    try:
        sys.argv["--in"] = snakemake.input[0]
        sys.argv["--count"] = snakemake.params.max_target_seqs
        sys.argv["--out"] = snakemake.output[0]
    except NameError as err:
        pass


def main():
    """Entry point."""
    try:
        parse_args()
        args = docopt(__doc__)
    except DocoptExit:
        raise DocoptExit
    try:
        lines = defaultdict(dict)
        chunk_counts = defaultdict(int)
        with open(args["--in"], "r") as fh:
            for line in fh.readlines():
                if "\n" not in line:
                    line += "\n"
                fields = line.split("\t")
                if fields[0]:
                    name, start = re.split("_-_", fields[0])
                    fields[0] = name
                    fields[3] = name
                    fields[9] = str(int(fields[9]) + int(start))
                    fields[10] = str(int(fields[10]) + int(start))
                    if start not in lines[name]:
                        lines[name][start] = []
                        chunk_counts[name] += 1
                    lines[name][start].append("\t".join(fields))
        with open(args["--out"], "w") as ofh:
            for name in chunk_counts.keys():
                length = len(lines[name])
                n = int(args["--count"])
                for start in lines[name].keys():
                    for i in range(n):
                        if i < len(lines[name][start]):
                            ofh.write("%s" % (lines[name][start][i]))
    except Exception as err:
        logger.error(err)
        exit(1)
