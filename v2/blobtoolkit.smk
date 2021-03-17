"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to run BlobToolKit sub-pipelines
-----------------------------------------

Requirements:
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p \
    --directory ~/workdir \
    --configfile example.yaml \
    -s blobtoolkit.smk
    -j 60

© 2021 Richard Challis (Wellcome Sanger Institute), MIT License
"""

import os

include: 'scripts/functions.py'

working_dir = os.getcwd()
parent_dir = os.path.dirname(working_dir)

rule all:
    """
    Dummy rule to define output
    """
    input:
        "%s/busco.stats" % parent_dir,
        "%s/minimap.stats" % parent_dir,
        "%s/diamond.stats" % parent_dir,
        "%s/blobtools.stats" % parent_dir,
        "%s/view.stats" % parent_dir

include: "rules/run_sub_pipeline.smk"
