#!/usr/bin/env python3

import os
import re
import yaml
import subprocess
import shlex
from pathlib import Path

COVERAGE = snakemake.input.cov
ASSEMBLY = snakemake.wildcards.assembly
OUTFILE = str(snakemake.output)
PORT = str(snakemake.params.port)
CLI = 'cli'
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
    cmds.append("%s %s --port %s %s --out %s/" % (CLI, ASSEMBLY, PORT, view, ASSEMBLY))

cmds.append("%s add --key static_plots=true ./%s" % (BLOBTOOLS, ASSEMBLY))

for cmd in cmds:
    subprocess.run(shlex.split(cmd), encoding='utf-8')

for filename in os.listdir(ASSEMBLY):
    p = Path("%s/%s" % (ASSEMBLY,filename))
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
