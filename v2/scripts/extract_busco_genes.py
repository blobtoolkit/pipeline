#!/usr/bin/env python3
"""Extract BUSCO gene sequences and set headers."""

import codecs
import logging
import re
import tarfile
from pathlib import Path

from docopt import docopt
from tolkein import tofile

docs = """
Extract BUSCO genes.

Usage: ./extract_busco_genes.py [--busco PATH...] [--out FASTA]

Options:
    --busco PATH         BUSCO full summary tsv file.
    --out FASTA          output FASTA filename.
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
        args["--busco"] = snakemake.input.busco
        args["--out"] = snakemake.output.fasta
    except NameError as err:
        logger.info(err)
        logger.info("Parsing parameters from command line")
    return args


if __name__ == "__main__":
    args = parse_args(docopt(docs))
    file_pattern = re.compile(r"busco_sequences\/(\w+?)_\w+\/(\d+at\d+).faa")
    header_pattern = re.compile(r">\S+:\d+-\d+")
    utf8reader = codecs.getreader("utf-8")
    try:
        with open(args["--out"], "w") as ofh:
            for busco_file in args["--busco"]:
                busco_dir = Path(busco_file).absolute().parent
                busco_seqs = "%s/busco_sequences.tar.gz" % busco_dir
                tar = tarfile.open(busco_seqs)
                for tarinfo in tar.getmembers():
                    if tarinfo.name.endswith(".faa"):
                        match = file_pattern.match(tarinfo.name)
                        status, busco_id = match.groups()
                        with utf8reader(
                            tar.extractfile(
                                tarinfo,
                            )
                        ) as fh:
                            header = next(fh).strip()
                            title = ""
                            if header_pattern.match(header):
                                title = "%s=%s=%s" % (header, busco_id, status)
                            else:
                                parts = header.split(" # ")
                                title = "%s:%s-%s=%s=%s" % (
                                    re.sub(r"_[^_]+$", "", parts[0]),
                                    parts[1],
                                    parts[2],
                                    busco_id,
                                    status,
                                )
                            ofh.write("%s\n" % title)
                            ofh.writelines(fh)
    except Exception as err:
        logger.error(err)
        exit(1)
