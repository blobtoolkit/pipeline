import os

"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to run BUSCO
--------------------------------------------------

Requirements:
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p \
    --directory ~/workdir \
    --configfile example.yaml \
    -s busco.smk
    -j 8

© 2021 Richard Challis (Wellcome Sanger Institute), MIT License
"""

include: 'scripts/functions.py'

rule all:
    """
    Dummy rule to define all outputs
    """
    input:
        expand("%s.busco.{lineage}.tsv" % config["assembly"]["prefix"], lineage=config['busco']['lineages'])


include: 'rules/run_busco5.smk'
include: 'rules/unzip_assembly_fasta.smk'
