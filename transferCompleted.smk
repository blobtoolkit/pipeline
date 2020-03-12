"""
https://github.com/blobtoolkit/insdc-pipeline

Pipeline to validate and transfer BlobDir datasets
----------------------------------------------

Requirements:
 - BlobTools2 (https://github.com/blobtoolkit/blobtools2)
 - Specification (https://github.com/blobtoolkit/specification)
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p --use-conda
    --directory path/to/workdir/
    --configfile path/to/config.yaml \
    -s transferCompleted.smk
    -j 8

© 2018-19 Richard Challis (University of Edinburgh), MIT License
© 2019-20 Richard Challis (Wellcome Sanger Institute), MIT License
"""

singularity: "docker://genomehubs/blobtoolkit:1.1"

include: 'scripts/functions.py'

supported_version = check_version()

use_singularity = check_config()

reads = get_read_info(config)
keep = False
if 'keep_intermediates' in config:
    keep = bool(config['keep_intermediates'])
asm = config['assembly']['prefix']
destdir = config['destdir']
rev = ''
if 'revision' in config:
    if config['revision'] > 0:
        rev = '.'+str(config['revision'])

rule all:
    """
    Dummy rule to set target of pipeline
    """
    input:
        "%s%s.complete" % (asm, rev)


include: 'rules/validate_dataset.smk'
include: 'rules/generate_images.smk'
include: 'rules/generate_summary.smk'
include: 'rules/checksum_files.smk'
include: 'rules/transfer_dataset.smk'
