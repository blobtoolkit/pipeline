#!/usr/bin/env python3

import os
import re
import yaml
import subprocess
from functions import get_read_info

CONFIG = snakemake.input.yaml
ASSEMBLY = snakemake.wildcards.assembly
INFILES = snakemake.input.stats
OUTFILE = str(snakemake.output)
GITDIR = snakemake.params.gitdir

meta = {}
with open(CONFIG,'r') as fh:
    config = yaml.load(fh)
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
