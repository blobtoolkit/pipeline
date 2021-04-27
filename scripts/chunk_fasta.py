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
from tolkein import tofile

docs = """
Chunk FASTA.

Usage: ./chunk_fasta.py [--in FASTA] [--chunk INT] [--overlap INT] [--max-chunks INT]
                        [--busco TSV] [--min-length INT] [--out FASTA] [--bed BEDFILE]

Options:
    --in FASTA           input FASTA file.
    --busco TSV          BUSCO full summary tsv file.
    --chunk INT          sequences greater than CHUNK bp will be split. [Default: 100000]
    --overlap INT        length of overlap when splitting sequences. [Default: 0]
    --max-chunks INT     maximum number of chunks to split a sequence into. [Default: 10]
    --min-length INT     minimum sequence length. [Default: 1000]
    --out FASTA          output FASTA filename or suffix. [Default: .chunked]
    --bed BEDFILE        output BED filename or suffix. [Default: .bed]
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


def chunk_fasta(
    fastafile, *, chunk=math.inf, overlap=0, max_chunks=math.inf, min_length=1000
):
    """Read FASTA file one sequence at a time and split long sequences into chunks."""
    cmd = "cat %s" % fastafile
    if fastafile.endswith(".gz"):
        cmd = "zcat %s" % fastafile
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
            elif seq_length > min_length:
                yield {
                    "title": title,
                    "seq": seq,
                    "chunks": 1,
                    "start": 0,
                    "end": seq_length,
                    "length": seq_length,
                }


def check_for_unmasked_bases(seq, min_unmasked=20):
    """Check sequence has runs of unmasked bases."""
    return bool(re.search(r"[ACGT]{" + str(min_unmasked) + "}", seq))


def check_for_masked_bases(seq, max_masked=20):
    """Check sequence has runs of masked bases."""
    return bool(re.search(r"[acgtnN]{" + str(max_masked) + "}", seq))


def parse_busco_full_summary(busco_file, chunk=100000):
    """Parse a busco full summary file."""
    logger.info("Parsing BUSCO full summary file")
    locations = defaultdict(list)
    with tofile.open_file_handle(busco_file) as fh:
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
    # with open("/tmp/data.tsv", "w") as fh:
    #     for tup in windows["LR862382.1"]:
    #         fh.write("\t".join([str(tup[0]), str(tup[2])]) + "\n")
    # quit()
    return windows


def chunk_by_busco(seq, args):
    """Use BUSCO gene density to define chunks."""
    subseq_found = False
    has_busco = False
    if seq["title"] in busco_windows:
        for window in busco_windows[seq["title"]]:
            if seq["start"] <= window[0] and seq["end"] >= window[1]:
                has_busco = True
                subseq = seq["seq"][window[0] - seq["start"] : window[1] - seq["start"]]
                # check there are runs of unmasked bases
                if check_for_masked_bases(subseq, 1000):
                    seq["seq"] = subseq
                    seq["start"] = window[0]
                    subseq_found = True
                    break
    if not has_busco:
        # extract subseq from middle of chunk
        chunk = int(args["--chunk"])
        offset = int((seq["length"] - chunk) / 2)
        if offset > 0:
            midseq = seq["seq"][offset:-offset]
            # check there are runs of unmasked bases
            if check_for_masked_bases(midseq, 1000):
                seq["seq"] = midseq
                seq["start"] += offset
                seq["end"] -= offset
                subseq_found = True
            else:
                # walk along sequence to find a chunk with fewer masked bases
                while offset > chunk:
                    offset -= chunk
                    lseq = seq["seq"][offset : offset + chunk]
                    if check_for_masked_bases(lseq, 1000):
                        seq["seq"] = lseq
                        seq["start"] += offset
                        seq["end"] = seq["start"] + chunk
                        subseq_found = True
                        break
                    rseq = seq["seq"][-offset - chunk : -offset]
                    if check_for_masked_bases(rseq, 1000):
                        seq["seq"] = rseq
                        seq["start"] = seq["end"] - offset - chunk
                        seq["end"] -= offset
                        subseq_found = True
                        break
    if subseq_found:
        seqs.append((seq))


def seq_stats(seq):
    """Calculate basic sequence statistics."""
    counter = Counter(seq)
    gc = 0
    at = 0
    n = 0
    masked = 0
    for base in ["g", "c", "G", "C"]:
        if base in counter:
            gc += counter[base]
    for base in ["a", "t", "A", "T"]:
        if base in counter:
            at += counter[base]
    for base in ["n", "N"]:
        if base in counter:
            n += counter[base]
    for base in ["a", "c", "g", "t", "n"]:
        if base in counter:
            masked += counter[base]
    atgc = at + gc
    if atgc > 0:
        gc = gc / atgc
    else:
        gc = 0
    return {
        "gc": gc,
        "n": n / len(seq),
        "ncount": n,
        "masked": masked / len(seq),
        "length": len(seq),
    }


def write_bedfiles(bed_data, args):
    """Write bedfiles."""
    if "--bed" not in args or not args["--bed"]:
        return
    lines = defaultdict(list)
    for title, arr in bed_data.items():
        sums = defaultdict(float)
        positions = []
        start = arr[0]["start"]
        end = arr[0]["end"]
        for obj in arr:
            location = "%s\t%d\t%d" % (title, obj["start"], obj["end"])
            for key, value in obj["stats"].items():
                start = min(obj["start"], start)
                end = max(obj["end"], end)
                if key != "ncount" and key != "length":
                    sums[key] += value * (obj["end"] - obj["start"])
                    lines["%s_windows" % key].append(
                        "%s\t.\t%.4f\n" % (location, value)
                    )
                else:
                    sums[key] += value
                    if key == "length":
                        position = sums[key]
                        lines["length_windows"].append(
                            "%s\t.\t%d\n" % (location, value)
                        )
                        lines["position_windows"].append(
                            "%s\t.\t%d\n" % (location, position)
                        )
                        positions.append(position)
                    else:
                        lines["%s_windows" % key].append(
                            "%s\t.\t%d\n" % (location, value)
                        )
                if key == "gc":
                    lines["mask"].append("%s\twindow\n" % location)
        location = "%s\t%d\t%d" % (title, start, end)
        length = end - start
        lines["mask"].append("%s\tfull\n" % location)
        for key, value in sums.items():
            if key != "ncount" and key != "length":
                lines[key].append("%s\t.\t%.4f\n" % (location, value / length))
            elif key == "length":
                lines["length"].append("%s\t.\t%d\n" % (location, value))
                lines["position"].append("%s\t.\t%d\n" % (location, value))
                lines["proportion"].append("%s\t.\t%.4f\n" % (location, value / length))
                for position in positions:
                    lines["proportion_windows"].append(
                        "%s\t.\t%.4f\n" % (location, position / length)
                    )
            else:
                lines[key].append("%s\t.\t%.d\n" % (location, value))
    for key, rows in lines.items():
        bedfile = args["--bed"]
        bedfile = re.sub(r"\.bed$", "", bedfile)
        bedfile += ".%s.bed" % key
        if bedfile.startswith("."):
            bedfile = "%s%s" % (args["--in"], bedfile)
        with open(bedfile, "w") as ofh:
            ofh.writelines(rows)


def parse_args(args):
    """Parse snakemake args if available."""
    try:
        args["--in"] = snakemake.input.fasta
        args["--chunk"] = int(snakemake.params.chunk)
        args["--overlap"] = int(snakemake.params.overlap)
        args["--max-chunks"] = int(snakemake.params.max_chunks)
        args["--min-length"] = int(snakemake.params.min_length)
        try:
            args["--busco"] = snakemake.input.busco
        except AttributeError:
            args["--busco"] = None
        try:
            args["--out"] = snakemake.output.fasta
        except AttributeError:
            args["--out"] = None
        try:
            args["--bed"] = snakemake.params.bed
        except AttributeError:
            args["--bed"] = None
    except NameError as err:
        logger.info(err)
        logger.info("Parsing parameters from command line")
    return args


if __name__ == "__main__":
    try:
        args = parse_args(docopt(docs))
        busco_windows = {}
        if "--busco" in args and args["--busco"] is not None:
            busco_windows = parse_busco_full_summary(args["--busco"])
        logger.info("Splitting %s into chunks" % args["--in"])
        bed_data = defaultdict(list)
        seqs = []
        for seq in chunk_fasta(
            args["--in"],
            chunk=int(args["--chunk"]),
            overlap=int(args["--overlap"]),
            max_chunks=int(args["--max-chunks"]),
            min_length=int(args["--min-length"]),
        ):
            if (
                busco_windows
                and seq["chunks"] == int(args["--max-chunks"])
                and seq["length"] > int(args["--chunk"]) + int(args["--overlap"])
            ):
                chunk_by_busco(seq, args)
            else:
                if "--bed" in args and args["--bed"]:
                    stats = seq_stats(seq["seq"])
                    bed_data[seq["title"]].append(
                        {
                            "start": seq["start"],
                            "end": seq["start"] + len(seq["seq"]),
                            "stats": stats,
                        }
                    )
                if check_for_unmasked_bases(seq["seq"]):
                    seqs.append((seq))
        if bed_data:
            write_bedfiles(bed_data, args)
        if args["--out"] is not None:
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
