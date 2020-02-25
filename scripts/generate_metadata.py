#!/usr/bin/env python3

import os
import re
import yaml
import subprocess
import logging
from functions import get_read_info

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
    # with open(CONFIG,'r') as fh:
    config = yaml.load(CONFIG)
    meta['assembly'] = config['assembly']
    meta['taxon'] = config['taxon']
    meta['settings'] = config['settings']
    meta['similarity'] = config['similarity']
    meta['reads'] = get_read_info(config)


    p = subprocess.Popen(['git','--git-dir',GITDIR,'rev-parse','--short','HEAD'], stdout=subprocess.PIPE,encoding='utf-8')
    (output, err) = p.communicate()
    p_status = p.wait()
    meta['settings']['commit'] = output.rstrip('\n')
    p = subprocess.Popen(['git','--git-dir',GITDIR,'remote','-v'], stdout=subprocess.PIPE,encoding='utf-8')
    (output, err) = p.communicate()
    p_status = p.wait()
    meta['settings']['pipeline'] = output.split()[1]

    for file in INFILES:
        accession = re.search(r'\.(.+).bam.stats',file.split(ASSEMBLY)[1]).group(1)
        sections = ('Total reads','Mapped reads','Proper-pairs','Both pairs mapped',
            'Average insert size', 'Median insert size')
        with open(file,'r') as fh:
            for line in fh:
                line = line.rstrip('\n')
                for section in sections:
                    if line.startswith(section):
                        meta['reads'][accession][section] = int(float(re.search(r'([\d\.]+)',line).group(1)))

    with open(OUTFILE,'w') as fh:
        fh.write(yaml.dump(meta,indent=1))

except Exception as err:
    logger.error(err)
    exit(1)
