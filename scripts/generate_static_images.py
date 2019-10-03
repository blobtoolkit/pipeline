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
             "--view snail --format png --format svg",
             "--stats"]

    if not COVERAGE:
        views = views[3:]

    cmds = []

    for view in views:
        cmds.append("%s %s --host %s %s --out %s/" % (CLI, ASSEMBLY, PORT, view, ASSEMBLY))

    cmds.append("%s add --key static_plots=true ./%s" % (BLOBTOOLS, ASSEMBLY))

    for cmd in cmds:
        subprocess.run(shlex.split(cmd), encoding='utf-8')

    for filename in os.listdir(ASSEMBLY):
        p = Path("%s/%s" % (ASSEMBLY, filename))
        parts = filename.split('.')
        if len(parts) == 3 and parts[0] == ASSEMBLY:
            if parts[2] == 'png' or parts[2] == 'svg':
                new_p = Path("%s/%s.%s" % (p.parent.as_posix(), parts[1], parts[2]))
                p.rename(new_p)
            elif parts[1] == 'current' and parts[2] == 'json':
                new_p = Path("%s/summary.%s" % (p.parent.as_posix(), parts[2]))
                p.rename(new_p)
        elif len(parts) == 4 and parts[0] == ASSEMBLY:
            if parts[3] == 'png' or parts[3] == 'svg':
                new_p = Path("%s/%s.%s.%s" % (p.parent.as_posix(), parts[1], parts[2], parts[3]))
                p.rename(new_p)
except Exception as err:
    logger.error(err)
    exit(1)
