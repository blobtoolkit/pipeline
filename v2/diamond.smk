import os

"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to run Diamond blastx
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

include: "scripts/functions.py"

busco_path = "../busco"

rule all:
    """
    Dummy rule to define output
    """
    input:
        "%s.diamond.reference_proteomes.out" % config["assembly"]["prefix"]


include: "rules/chunk_fasta_by_busco.smk"
include: "rules/run_diamond_blastx.smk"
include: "rules/unchunk_blast.smk"
