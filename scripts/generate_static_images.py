#!/usr/bin/env python3

import glob
import logging
import os
import shlex
import signal
import subprocess
from pathlib import Path

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
    TIMEOUT = int(snakemake.params.timeout)
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
            "%s view --host %s --timeout %d --ports %s %s --out ./%s/ %s"
            % (BLOBTOOLS, HOST, TIMEOUT, PORTS, view, BLOBDIR, BLOBDIR)
        )

    cmds.append(
        "%s filter --summary %s.summary.json ./%s" % (BLOBTOOLS, BLOBDIR, BLOBDIR)
    )

    cmds.append("%s add --key static_plots=true ./%s" % (BLOBTOOLS, BLOBDIR))

    for cmd in cmds:
        logger.info(cmd)
        with subprocess.Popen(
            shlex.split(cmd),
            stdout=subprocess.PIPE,
            preexec_fn=os.setsid,
            encoding="utf-8",
        ) as process:
            try:
                output = process.communicate(timeout=1800)[0]
            except subprocess.TimeoutExpired:
                os.killpg(
                    process.pid, signal.SIGINT
                )  # send signal to the process group
                output = process.communicate()[0]
        # try:
        #     subprocess.check_output(shlex.split(cmd), encoding="utf-8", timeout=1800)
        # except Exception as err:
        #     raise err

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
    for pngpath in glob.iglob(os.path.join(BLOBDIR, "%s.*.png" % BLOBDIR)):
        os.remove(pngpath)
    for svgpath in glob.iglob(os.path.join(BLOBDIR, "%s.*.svg" % BLOBDIR)):
        os.remove(svgpath)
    for jsonpath in glob.iglob(os.path.join(BLOBDIR, "%s.*.json" % BLOBDIR)):
        os.remove(jsonpath)
    exit(1)
