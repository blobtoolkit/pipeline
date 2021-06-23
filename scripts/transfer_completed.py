#!/usr/bin/env python3

"""
Transfer completed BlobDirs and intermediate files.

Usage:
  transfer_completed.py [--in PATH] [--out PATH] [--bin PATH]

Options:
  --in PATH              Path to input directory [default: .]
  --out PATH             Path to output directory [default: ./complete]
  --bin PATH             Path to bin directory. If specified, unwanted files will
                         be moved here rather than deleted.
"""

import glob
import gzip
import shutil
import tarfile
from pathlib import Path

import yaml
from docopt import docopt
from tolkein import tofetch, tofile, tolog

LOGGER = tolog.logger(__name__)


def compress_file(infile, outfile):
    """Compress a file using gzip."""
    LOGGER.info("Compressing %s to %s", infile, outfile)
    with gzip.open(outfile, "wb") as ofh:
        with open(infile, "rb") as ifh:
            shutil.copyfileobj(ifh, ofh)


def compress_content(indir, destdir):
    """Compress all files in directory."""
    Path(destdir).mkdir(parents=True, exist_ok=True)
    for infile in glob.glob(r"%s/*" % indir):
        outfile = "%s/%s.gz" % (destdir, Path(infile).name)
        compress_file(infile, outfile)


def tar_directory(indir, destdir, *, compress=False):
    LOGGER.info("Archiving %s to %s.tar", indir, destdir)
    if compress:
        compress_content(indir, destdir)
    else:
        shutil.copytree(indir, destdir)
    with tarfile.open(destdir + ".tar", mode="w") as archive:
        archive.add(destdir, Path(destdir).name, recursive=True)
    shutil.rmtree(destdir)


def transfer_files(pattern, destdir, *, compress=False, rename=False):
    """Transfer files to output directory."""
    for file in glob.glob(pattern):
        filepath = Path(file)
        destfile = "%s/%s" % (destdir, filepath.name)
        if rename:
            destfile = destfile.replace(rename[0], rename[1])
        if filepath.is_dir():
            tar_directory(file, destfile, compress=compress)
            continue
        if compress:
            destfile += ".gz"
            compress_file(file, destfile)
            continue
        LOGGER.info("Moving %s to %s", file, destdir)
        shutil.move(file, destdir)


def create_static_directory(indir):
    """Move image files to static directory."""
    staticdir = "%s_static" % indir
    Path(staticdir).mkdir(parents=True, exist_ok=True)
    for file in glob.glob(r"%s/*.png" % indir):
        shutil.move(file, staticdir)


def remove_unwanted_files(indir, bindir):
    """Move unwanted files to bin or delete."""
    LOGGER.info("Removing unwanted files")
    if bindir is not None:
        LOGGER.info("Moving %s to %s", indir, bindir)
        Path(bindir).mkdir(parents=True, exist_ok=True)
        shutil.move(indir, bindir)
    else:
        LOGGER.info("Deleting %s")
        shutil.rmtree(indir)


if __name__ == "__main__":
    opts = docopt(__doc__)
    data = tofile.read_file("%s/config.yaml" % opts["--in"])
    meta = yaml.full_load(data)
    prefix = meta["assembly"]["prefix"]
    indir = opts["--in"]
    outdir = "%s/%s" % (opts["--out"], prefix)
    bindir = opts.get("--bin", None)
    Path(outdir).mkdir(parents=True, exist_ok=True)
    create_static_directory("%s/view/%s" % (indir, prefix))
    transfer_files(r"%s/blastn/%s.*.out" % (indir, prefix), outdir, compress=True)
    transfer_files("%s/blobtools/%s.meta.yaml" % (indir, prefix), outdir, compress=True)
    transfer_files(r"%s/busco/%s.busco.*_odb10" % (indir, prefix), outdir)
    transfer_files("%s/config.yaml" % indir, outdir)
    transfer_files(
        "%s/cov_stats/%s.chunk_stats.tsv" % (indir, prefix), outdir, compress=True
    )
    transfer_files(
        r"%s/diamond/%s.*.out" % (indir, prefix),
        outdir,
        rename=("diamond.", "diamond_blastx."),
        compress=True,
    )
    transfer_files(
        "%s/diamond/%s.fasta.chunks" % (indir, prefix),
        outdir,
        rename=("fasta.chunks", "chunks.fasta"),
        compress=True,
    )
    transfer_files(
        r"%s/diamond_blastp/%s.*.out" % (indir, prefix),
        outdir,
        rename=("diamond.", "diamond_blastp."),
        compress=True,
    )
    transfer_files(r"%s/minimap/%s.*.bam" % (indir, prefix), outdir)
    transfer_files("%s/view/%s" % (indir, prefix), outdir, compress=True)
    transfer_files("%s/view/%s_static" % (indir, prefix), outdir)
    transfer_files(
        "%s/window_stats/%s.window_stats.tsv" % (indir, prefix), outdir, compress=True
    )
    transfer_files(
        "%s/windowmasker/%s.windowmasker.fasta" % (indir, prefix), outdir, compress=True
    )
    remove_unwanted_files(indir, bindir)
