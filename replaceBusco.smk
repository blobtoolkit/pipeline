import os

"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to run replace BUSCO results in a BlobDir
--------------------------------------------------

Requirements:
 - BlobTools2 (https://github.com/blobtoolkit/blobtools2)
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p --use-conda \
    --directory ~/workdir \
    --configfile example.yaml \
    -s replaceHits.smk
    -j 8

© 2018-19 Richard Challis (University of Edinburgh), MIT License
© 2019-20 Richard Challis (Wellcome Sanger Institute), MIT License
"""

singularity: "docker://genomehubs/blobtoolkit:1.1"

include: 'scripts/functions.py'

use_singularity = check_config()

supported_version = check_version()

multicore = int(os.getenv('MULTICORE', 16))
maxcore = int(os.getenv('MAXCORE', 32))

similarity = apply_similarity_search_defaults()

reads = get_read_info(config)
keep = False
keep = False
if 'keep_intermediates' in config:
    keep = bool(config['keep_intermediates'])
asm = config['assembly']['prefix']
rev = ''
if 'revision' in config:
    if config['revision'] > 0:
        rev = '.'+str(config['revision'])

rule all:
    """
    Dummy rule to set blobDB as target of pipeline
    """
    input:
        "%s.busco.removed" % asm,
        expand("%s/{lineage}_busco.json" % asm,lineage=config['busco']['lineages'])


include: 'rules/fetch_busco_lineage.smk'
include: 'rules/run_busco.smk'
include: 'rules/blobtoolkit_remove_busco.smk'
include: 'rules/blobtoolkit_add_busco.smk'
