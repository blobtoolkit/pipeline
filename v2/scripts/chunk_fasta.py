#!/usr/bin/env python3
"""Split long sequences into chunks."""

import glob
import logging
import math
import os.path
import re
import shlex
import sys
import time
from collections import Counter, defaultdict
from itertools import groupby
from random import shuffle
from subprocess import PIPE, Popen

from docopt import docopt

docs = """
Chunk FASTA.

Usage: ./chunk_fasta.py [--in FASTA] [--chunk INT] [--overlap INT] [--max-chunks INT]
                        [--busco TSV] [--out CHUNKED_FASTA]

Options:
    --in FASTA           input FASTA file.
    --busco TSV          BUSCO full summary tsv file.
    --chunk INT          sequences greater than CHUNK bp will be split. [Default: 100000]
    --overlap INT        length of overlap when splitting sequences. [Default: 500]
    --max-chunks INT     maximum number of chunks to split a sequence into. [Default: 10]
    --out CHUNKED_FASTA  output filename. [Default: .chunked]
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


def chunk_size(value):
    """Calculate nice value for chunk size."""
    mag = math.floor(math.log10(value))
    first = int(str(value)[:2]) + 1
    chunk_size = first * pow(10, mag - 1)
    return chunk_size


def chunk_fasta(fastafile, chunk=math.inf, overlap=0, max_chunks=math.inf):
    """Read FASTA file one sequence at a time and split long sequences into chunks."""
    cmd = "cat %s" % fastafile
    # TODO: read gzipped files if needed
    # cmd = "pigz -dc %s" % fastafile
    title = ""
    seq = ""
    segment = chunk + overlap
    with Popen(shlex.split(cmd), encoding="utf-8", stdout=PIPE, bufsize=4096) as proc:
        faiter = (x[1] for x in groupby(proc.stdout, lambda line: line[0] == ">"))
        for header in faiter:
            title = header.__next__()[1:].strip().split()[0]
            seq = "".join(map(lambda s: s.strip(), faiter.__next__()))
            seq_length = len(seq)
            my_chunk = chunk
            if seq_length > segment:
                n = (seq_length + chunk) // chunk
                if n > max_chunks:
                    my_chunk = chunk_size(seq_length / max_chunks)
                    n = max_chunks
                for i in range(0, seq_length, my_chunk):
                    subseq = seq[i : i + my_chunk + overlap]
                    yield {
                        "title": title,
                        "seq": subseq,
                        "chunks": n,
                        "start": i,
                        "end": i + my_chunk + overlap,
                        "length": my_chunk + overlap,
                    }
            else:
                yield {
                    "title": title,
                    "seq": seq,
                    "chunks": 1,
                    "start": 0,
                    "end": seq_length,
                    "length": seq_length,
                }


def check_for_unmasked_bases(seq, min_unmasked=12):
    """Check sequence has a minimum number of unmasked bases."""
    return bool(re.search(r"[ACGT]{20}", seq))


def parse_busco_full_summary(busco_file, chunk=100000):
    """Parse a busco full summary file."""
    logger.info("Parsing BUSCO full summary file")
    locations = defaultdict(list)
    with open(busco_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[1] in ("Complete", "Duplicated"):
                locations[parts[2]].append((int(parts[3]), int(parts[4])))
    windows = {}
    for title, tuples in locations.items():
        tuples.sort(key=lambda tup: tup[0])
        windows[title] = []
        start_index = 0
        for location in tuples:
            windows[title].append([location[0], location[0] + chunk, 0])
            for window in windows[title][start_index:]:
                if location[1] < window[1]:
                    window[2] += 1
                else:
                    start_index += 1
        windows[title].sort(key=lambda window: window[2], reverse=True)
    return windows


if __name__ == "__main__":
    args = docopt(docs)
    try:
        args["--in"] = snakemake.input.fasta
        args["--busco"] = snakemake.input.busco
        args["--chunk"] = int(snakemake.params.chunk)
        args["--overlap"] = int(snakemake.params.overlap)
        args["--max-chunks"] = int(snakemake.params.max_chunks)
        args["--out"] = snakemake.output[0]
    except NameError as err:
        logger.info(err)
        logger.info("Parsing parameters from command line")
    busco_windows = {}
    if "--busco" in args:
        busco_windows = parse_busco_full_summary(args["--busco"])
    logger.info("Splitting %s into chunks" % args["--in"])
    try:
        seqs = []
        for seq in chunk_fasta(
            args["--in"],
            chunk=int(args["--chunk"]),
            overlap=int(args["--overlap"]),
            max_chunks=int(args["--max-chunks"]),
        ):
            if (
                busco_windows
                and seq["chunks"] == int(args["--max-chunks"])
                and seq["length"] > int(args["--chunk"]) + int(args["--overlap"])
            ):
                has_busco = False
                if seq["title"] in busco_windows:
                    for window in busco_windows[seq["title"]]:
                        if seq["start"] <= window[0] and seq["end"] >= window[1]:
                            has_busco = True
                            seq["seq"] = seq["seq"][
                                window[0] - seq["start"] : window[1] - seq["start"]
                            ]
                            seq["start"] = window[0]
                            break
                if not has_busco:
                    offset = int((seq["length"] - int(args["--chunk"])) / 2)
                    if offset > 0:
                        seq["seq"] = seq["seq"][offset:-offset]
                        seq["start"] += offset
                        seq["end"] -= offset
            if check_for_unmasked_bases(seq["seq"]):
                seqs.append((seq))
        chunked = ""
        for seq in seqs:
            chunked += ">%s_-_%d\n" % (seq["title"], seq["start"])
            chunked += "%s\n" % seq["seq"]
        outfile = args["--out"]
        if outfile.startswith("."):
            outfile = "%s%s" % (args["--in"], args["--out"])
        with open(outfile, "w") as ofh:
            ofh.writelines(chunked)
    except Exception as err:
        logger.error(err)
        exit(1)