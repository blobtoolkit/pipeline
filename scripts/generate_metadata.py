#!/usr/bin/env python3

import os
import re
import yaml
import subprocess
import logging
from collections import OrderedDict

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
    CONFIG = snakemake.config
    ASSEMBLY = snakemake.wildcards.assembly
    INFILES = snakemake.input.stats
    OUTFILE = str(snakemake.output)
    GITDIR = snakemake.params.gitdir
except Exception as err:
    logger.error(err)
    exit(1)

try:
    meta = {}
    config = CONFIG
    meta['assembly'] = config['assembly']
    meta['taxon'] = config['taxon']
    meta['settings'] = config['settings']
    meta['similarity'] = config['similarity']

    meta['reads'] = {}
    strategies = ['paired', 'single']
    for strategy in strategies:
        for arr in config['reads'][strategy]:
            if arr:
                meta['reads'][arr[0]] = {
                    'strategy': strategy,
                    'platform': arr[1],
                }
                if len(arr) == 4:
                    meta['reads'][arr[0]]['url'] = re.split(',|;', arr[3])

    p = subprocess.Popen(['git', '--git-dir', GITDIR, 'rev-parse', '--short', 'HEAD'],
                         stdout=subprocess.PIPE, encoding='utf-8')
    (output, err) = p.communicate()
    p_status = p.wait()
    meta['settings']['commit'] = output.rstrip('\n')
    p = subprocess.Popen(['git', '--git-dir', GITDIR, 'remote', '-v'],
                         stdout=subprocess.PIPE,
                         encoding='utf-8')
    (output, err) = p.communicate()
    p_status = p.wait()
    meta['settings']['pipeline'] = output.split()[1]

    for file in INFILES:
        accession = re.search(r'\.(.+).bam.stats', file.split(ASSEMBLY)[1]).group(1)
        sections = ('Total reads', 'Mapped reads', 'Proper-pairs', 'Both pairs mapped',
                    'Average insert size', 'Median insert size')
        with open(file, 'r') as fh:
            for line in fh:
                line = line.rstrip('\n')
                for section in sections:
                    if line.startswith(section):
                        meta['reads'][accession][section] = int(float(re.search(r'([\d\.]+)', line).group(1)))

    yaml.add_representer(OrderedDict,
                         lambda dumper,
                         data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
    yaml.add_representer(tuple,
                         lambda dumper,
                         data: dumper.represent_sequence('tag:yaml.org,2002:seq', data))
    yaml.Dumper.ignore_aliases = lambda *args: True
    with open(OUTFILE, 'w') as fh:
        fh.write(yaml.dump(meta, indent=1))

except Exception as err:
    logger.error(err)
    exit(1)
