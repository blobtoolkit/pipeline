#!/usr/bin/env python3
"""Count busco genes."""

import logging
import math
import statistics
from collections import defaultdict
from pathlib import Path

from docopt import docopt
from tolkein import tofile

docs = """
Count BUSCO genes.

Usage: ./count_busco_genes.py [--in TSV...] [--mask TSV] [--out TSV]

Options:
    --in TSV      chunked summary stats tsv file.
    --mask TSV    BED or BED-like TSV format mask file to specify sequence chunks.
    --out TSV     output TSV filename or suffix.
"""

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


def parse_args(args):
    """Parse snakemake args if available."""
    try:
        args["--in"] = snakemake.input.busco
        args["--mask"] = snakemake.input.mask
        args["--out"] = snakemake.output.tsv
    except NameError as err:
        logger.info(err)
        logger.info("Parsing parameters from command line")
    return args


def load_mask(filename):
    """Load bed file as mask."""
    mask = defaultdict(dict)
    header = []
    with tofile.open_file_handle(filename) as fh:
        for line in fh.readlines():
            seqid, start, end, *cols = line.rstrip().split("\t")
            if cols is None:
                cols = []
            if seqid == "sequence" and start == "start":
                header = cols
                continue
            mask[seqid].update({int(start): {"end": int(end), "cols": cols}})
    return mask, header


def parse_busco_summary(filename, mask, header):
    """Parse chunked values into dict."""
    with tofile.open_file_handle(filename) as fh:
        buscos = defaultdict(list)
        for line in fh.readlines():
            if line.startswith("#"):
                if line.startswith("# The lineage dataset is:"):
                    meta = line.split()
                    lineage = meta[5]
                    header.append("%s_count" % lineage)
                continue
            busco, status, *rest = line.rstrip().split("\t")
            if status in {"Fragmented", "Missing"}:
                continue
            seqid, start, *rest = rest
            buscos[seqid].append(int(start))
        for seqid in mask:
            starts = sorted(buscos[seqid])
            i = 0
            for start, obj in mask[seqid].items():
                ctr = 0
                while i < len(starts):
                    if starts[i] >= start:
                        if starts[i] > obj["end"]:
                            break
                        ctr += 1
                        i += 1
                obj["cols"].append(ctr)

            # if header is None:
            #     header = {key: idx + 3 for idx, key in enumerate(row[3:])}
            #     continue
            # seqid = row[0]
            # chunk_length = int(row[2]) - int(row[1])
            # if chunk_length > interval:
            #     interval = chunk_length
            # lengths[seqid] += chunk_length
            # for key, idx in header.items():
            #     values[seqid][key].append(float(row[idx]))
    return mask, header


if __name__ == "__main__":
    try:
        args = parse_args(docopt(docs))
        mask, header = load_mask(args["--mask"])
        for buscofile in args["--in"]:
            mask, header = parse_busco_summary(buscofile, mask, header)
        outfile = args["--out"]
        if outfile.startswith("."):
            outfile = "%s%s" % (args["--in"], args["--out"])
        header = ["sequence", "start", "end"] + header
        rows = ["%s\n" % "\t".join(header)]
        for seqid in mask:
            for start, obj in mask[seqid].items():
                row = [seqid, start, obj["end"]] + obj["cols"]
                rows.append("%s\n" % "\t".join([str(v) for v in row]))
        with open(outfile, "w") as ofh:
            ofh.writelines(rows)
    except Exception as err:
        raise err
        logger.error(err)
        exit(1)
