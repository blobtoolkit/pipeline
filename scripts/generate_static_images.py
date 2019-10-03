#!/usr/bin/env python3

import os
import re
import yaml
import subprocess
import shlex
import logging
from pathlib import Path

logger_config = {
    'level': logging.INFO,
    'format': '%(asctime)s [%(levelname)s] line %(lineno)d %(message)s',
    'filemode': 'w'
}
try:
    logger_config.update({'filename': snakemake.log[0]})
except NameError as err:
    pass
logging.basicConfig(**logger_config)
logger = logging.getLogger()

try:
    COVERAGE = snakemake.input.cov
    ASSEMBLY = snakemake.wildcards.assembly
    OUTFILE = str(snakemake.output)
    HOST = str(snakemake.params.host)
    CLI = 'blobtools view'
    BLOBTOOLS = 'blobtools'

    views = [
             "--view blob --param plotShape=circle --format png --format svg",
             "--view blob --param plotShape=hex --format png --format svg",
             "--view blob --param plotShape=square --format png --format svg",
             "--view cumulative --format png --format svg",
             "--view snail --format png --format svg"]

    if not COVERAGE:
        views = views[3:]

    cmds = []

    for view in views:
        cmds.append("%s %s --host %s %s --out %s/" % (CLI, ASSEMBLY, HOST, view, ASSEMBLY))

    cmds.append("%s add --key static_plots=true ./%s" % (BLOBTOOLS, ASSEMBLY))

    for cmd in cmds:
        subprocess.run(shlex.split(cmd), encoding='utf-8')

    for filename in os.listdir(ASSEMBLY):
        p = Path("%s/%s" % (ASSEMBLY, filename))
        parts = filename.split('.')
        if parts[0] == ASSEMBLY:
            if parts[-1] == 'png' or parts[-1] == 'svg':
                if parts[2].isdigit():
                    parts = parts[2:]
                else:
                    parts = parts[1:]
                new_p = Path("%s/%s" % (p.parent.as_posix(), '.'.join(parts)))
                p.rename(new_p)
except Exception as err:
    logger.error(err)
    exit(1)
