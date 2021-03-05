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
    -s diamond.smk
    -j 8

© 2021 Richard Challis (Wellcome Sanger Institute), MIT License
"""

include: 'scripts/functions.py'

rule all:
    """
    Dummy rule to define output
    """
    input:
        expand("%s.{sra}.bam" % config["assembly"]["prefix"], sra=reads_by_prefix(config).keys())


include: "rules/make_diamond_db.smk"
include: "rules/chunk_fasta.smk"
include: "rules/run_diamond_blastx.smk"
include: "rules/unchunk_blast.smk"
