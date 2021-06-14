#!/usr/bin/env python3
"""Split long sequences into chunks."""

import logging
import math
import statistics
from collections import defaultdict
from pathlib import Path

from docopt import docopt
from tolkein import tofile

docs = """
Window stats.

Usage: ./window_stats.py [--in TSV] [--window FLOAT...] [--min-window-length INT] [--min-window-count INT] [--out TSV]

Options:
    --in TSV                 chunked summary stats tsv file.
    --min-window-length INT  minimum length of a window. [Default: 100000]
    --min-window-count INT   minimum number of windows. [Default: 1]
    --window INT             window size or proportion. [Default: 1]
    --out TSV                output TSV filename or suffix.
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
        args["--in"] = snakemake.input.tsv
        args["--window"] = snakemake.params.window
        args["--out"] = snakemake.output.tsv
    except NameError as err:
        logger.info(err)
        logger.info("Parsing parameters from command line")
    return args


def parse_chunked_values(filename):
    """Parse chunked values into dict."""
    interval = 0
    header = None
    values = defaultdict(lambda: defaultdict(list))
    lengths = defaultdict(int)
    with tofile.open_file_handle(filename) as fh:
        for line in fh.readlines():
            row = line.rstrip().split("\t")
            if header is None:
                header = {key: idx + 3 for idx, key in enumerate(row[3:])}
                continue
            seqid = row[0]
            chunk_length = int(row[2]) - int(row[1])
            if chunk_length > interval:
                interval = chunk_length
            lengths[seqid] += chunk_length
            for key, idx in header.items():
                values[seqid][key].append(float(row[idx]))
    return lengths, values, interval


def round_to_interval(value, interval):
    return round(
        round(value / interval + 0.5) * interval, -int(math.floor(math.log10(interval)))
    )


def get_window_size(length, interval, window, min_length, min_count):
    """Get size of window based on options."""
    if window == 1:
        return length
    if window < 1:
        window_length = round_to_interval(length * window, interval)
        if window_length > min_length:
            return window_length
        return None
    if length / window >= min_count:
        return round_to_interval(window, interval)
    return None


def calculate_mean(arr, log):
    """Calculate mean and sd of arr values."""
    n = len(arr)
    if log:
        logged_arr = [math.log10(value) for value in arr if value > 0]
        mean = math.pow(10, statistics.mean(logged_arr))
        sd = math.pow(10, statistics.stdev(logged_arr))
        return mean, sd, n
    mean = statistics.mean(arr)
    sd = statistics.stdev(arr)
    return mean, sd, n


def calculate_window_stats(lengths, chunks, window, interval, args):
    """Calculate mean and sd in windows across each sequence."""
    values = defaultdict(dict)
    for seqid, length in lengths.items():
        window_size = get_window_size(
            length,
            interval,
            window,
            int(args["--min-window-length"]),
            int(args["--min-window-count"]),
        )
        if window_size is not None:
            for key, arr in chunks[seqid].items():
                start_i = 0
                start_pos = 0
                while start_pos < length:
                    end_pos = min(start_pos + window_size, length)
                    # window_length = start_pos - end_pos
                    mid_pos = start_pos + (end_pos - start_pos) / 2
                    end_i = round(end_pos / interval + 0.5) - 1
                    if window == 1:
                        proportion = "1"
                    else:
                        proportion = "%.3f" % (mid_pos / length)
                    if start_pos not in values[seqid]:
                        values[seqid][start_pos] = {
                            "end": str(end_pos),
                            "proportion": proportion,
                        }
                    if key.endswith("count"):
                        values[seqid][start_pos][key] = "%d" % sum(arr[start_i:end_i])
                        # values[seqid][start_pos][
                        #     key.replace("_count", "_cpm")
                        # ] = "%.3f" % (sum(arr[start_i:end_i]) / window_length * 1000000)
                    else:
                        mean, sd, n = calculate_mean(
                            arr[start_i : end_i + 1], key.endswith("_cov")
                        )
                        values[seqid][start_pos][key] = "%.3f" % mean
                        values[seqid][start_pos]["%s_sd" % key] = "%.3f" % sd
                        values[seqid][start_pos]["%s_n" % key] = "%d" % n
                    start_pos = end_pos
                    start_i = end_i + 1
    return values


if __name__ == "__main__":
    try:
        args = parse_args(docopt(docs))
        lengths, chunks, interval = parse_chunked_values(args["--in"])
        outfile = args["--out"]
        if outfile.startswith("."):
            outfile = "%s%s" % (args["--in"], args["--out"])
        outfile = Path(outfile)
        suffix = outfile.suffix
        filename = outfile.stem
        for window in args["--window"]:
            values = calculate_window_stats(
                lengths, chunks, float(window), interval, args
            )
            rows = []
            header = None
            for seqid, obj in values.items():
                for start_pos, entry in obj.items():
                    if header is None:
                        header = ["sequence", "start", "end"]
                        for key in entry.keys():
                            if key not in header:
                                header.append(key)
                        rows.append("\t".join(header) + "\n")
                    row = [seqid, str(start_pos)]
                    for key in header[2:]:
                        row.append(entry[key])
                    rows.append("\t".join(row) + "\n")
            filetag = ""
            if window != 1:
                filetag = ".%s" % str(window)
            with open("%s%s%s" % (filename, filetag, suffix), "w") as fh:
                fh.writelines(rows)

        # if args["--out"] is not None:
        #     chunked = ""
        #     for seq in seqs:
        #         chunked += ">%s_-_%d\n" % (seq["title"], seq["start"])
        #         chunked += "%s\n" % seq["seq"]
        #     outfile = args["--out"]
        #     if outfile.startswith("."):
        #         outfile = "%s%s" % (args["--in"], args["--out"])
        #     with open(outfile, "w") as ofh:
        #         ofh.writelines(chunked)
    except Exception as err:
        logger.error(err)
        exit(1)
