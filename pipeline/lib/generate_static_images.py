#!/usr/bin/env python3

"""
Generate static image files.

Usage:
  generate-static-images --blobdir STRING [--coverage] [--host STRING]
                         [--ports RANGE] [--timeout INT] [--path PATH]

Options:
  --blobdir STRING  BlobDir to generate static views for.
  --coverage        Flag to indicate if coverage data are available and blob views should be generated. [Default False].
  --host STRING     Host name. [Default: http://localhost]
  --path STRING     Path to directory containing BlobDir. [Default: .]
  --ports RANGE     Range of ports available for hosting BlobTools Viewer/API. [Default: 8000-8099] = str(snakemake.params.ports)
  --timeout INT     Maximum time (in seconds) to wait for images. [Default: 30]
"""

import glob
import logging
import os
import shlex
import signal
import subprocess
import sys
from pathlib import Path

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
        sys.argv["--blobdir"] = snakemake.wildcards.blobdir
        sys.argv["--coverage"] = snakemake.input.cov
        sys.argv["--host"] = str(snakemake.params.host)
        sys.argv["--ports"] = str(snakemake.params.ports)
        sys.argv["--timeout"] = int(snakemake.params.timeout)
    except NameError as err:
        logger.info(err)
        logger.info("Parsing parameters from command line")


def main():
    """Entry point."""
    args = docopt(__doc__)

    blob_path = "%s/%s" % (args["--path"], args["--blobdir"])

    try:
        views = [
            "--view blob --param plotShape=circle --param largeFonts=true --format png",
            "--view blob --param plotShape=hex --param largeFonts=true --format png",
            "--view blob --param plotShape=square --param largeFonts=true --format png",
            "--view blob --param plotShape=kite --param largeFonts=true --format png",
            "--view cumulative --param largeFonts=true --format png",
            "--view snail --param largeFonts=true --format png",
        ]

        if not args["--coverage"]:
            views = views[4:]

        cmds = []

        for view in views:
            cmds.append(
                "blobtools view --host %s --timeout %d --ports %s %s --out %s/ %s"
                % (
                    args["--host"],
                    args["--timeout"],
                    args["--ports"],
                    view,
                    blob_path,
                    args["--blobdir"],
                )
            )

        cmds.append(
            "blobtools filter --summary %s.summary.json %s/%s"
            % (args["--blobdir"], blob_path)
        )

        cmds.append("blobtools add --key static_plots=true %s/%s" % (blob_path))

        for cmd in cmds:
            logger.info(cmd)
            with subprocess.Popen(
                shlex.split(cmd),
                stdout=subprocess.PIPE,
                preexec_fn=os.setsid,
                encoding="utf-8",
            ) as process:
                try:
                    process.communicate(timeout=1800)[0]
                except subprocess.TimeoutExpired:
                    os.killpg(
                        process.pid, signal.SIGINT
                    )  # send signal to the process group
                    process.communicate()[0]

        for filename in os.listdir(blob_path):
            p = Path("%s/%s" % (blob_path, filename))
            parts = filename.split(".")
            if filename.startswith(blob_path):
                if (
                    filename.endswith("png")
                    or filename.endswith("svg")
                    or filename.endswith("json")
                ):
                    if parts[1].isdigit():
                        parts = parts[2:]
                    else:
                        parts = parts[1:]
                    new_p = Path(
                        "%s/%s"
                        % (
                            p.parent.as_posix(),
                            filename.replace("%s." % args["--blobdir"], ""),
                        )
                    )
                    p.rename(new_p)
    except Exception as err:
        logger.error(err)
        for pngpath in glob.iglob(
            os.path.join(blob_path, "%s.*.png" % args["--blobdir"])
        ):
            os.remove(pngpath)
        for svgpath in glob.iglob(
            os.path.join(blob_path, "%s.*.svg" % args["--blobdir"])
        ):
            os.remove(svgpath)
        for jsonpath in glob.iglob(
            os.path.join(blob_path, "%s.*.json" % args["--blobdir"])
        ):
            os.remove(jsonpath)
        exit(1)


if __name__ == "__main__":
    main()
