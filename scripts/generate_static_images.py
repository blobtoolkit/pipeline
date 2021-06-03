#!/usr/bin/env python3

import logging
import os
import re
import shlex
import subprocess
from pathlib import Path

import yaml

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

try:
    COVERAGE = snakemake.input.cov
    BLOBDIR = snakemake.wildcards.blobdir
    OUTFILE = str(snakemake.output)
    HOST = str(snakemake.params.host)
    PORTS = str(snakemake.params.ports)
    BLOBTOOLS = "blobtools"

    views = [
        "--view blob --param plotShape=circle --param largeFonts=true --format png",
        "--view blob --param plotShape=hex --param largeFonts=true --format png",
        "--view blob --param plotShape=square --param largeFonts=true --format png",
        "--view blob --param plotShape=kite --param largeFonts=true --format png",
        "--view cumulative --param largeFonts=true --format png",
        "--view snail --param largeFonts=true --format png",
    ]

    if not COVERAGE:
        views = views[4:]

    cmds = []

    for view in views:
        cmds.append(
            "%s view --host %s --ports %s %s --out ./%s/ %s"
            % (BLOBTOOLS, HOST, PORTS, view, BLOBDIR, BLOBDIR)
        )

    cmds.append(
        "%s filter --summary %s.summary.json ./%s" % (BLOBTOOLS, BLOBDIR, BLOBDIR)
    )

    cmds.append("%s add --key static_plots=true ./%s" % (BLOBTOOLS, BLOBDIR))

    for cmd in cmds:
        logger.info(cmd)
        try:
            subprocess.run(shlex.split(cmd), encoding="utf-8")
        except Exception as err:
            logger.error(err)

    for filename in os.listdir(BLOBDIR):
        p = Path("%s/%s" % (BLOBDIR, filename))
        parts = filename.split(".")
        if filename.startswith(BLOBDIR):
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
                    % (p.parent.as_posix(), filename.replace("%s." % BLOBDIR, ""))
                )
                p.rename(new_p)
except Exception as err:
    logger.error(err)
    exit(1)
