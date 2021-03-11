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

rule all:
    """
    Dummy rule to define output
    """
    input:
        "%s.%s.out" % (config["assembly"]["prefix"], diamond_db_name(config))


# include: "rules/make_taxid_list.smk"
# include: "rules/make_masked_list.smk"
# include: "rules/make_diamond_db.smk"
# include: "rules/chunk_fasta.smk"
# include: "rules/run_diamond_blastx.smk"
# include: "rules/unchunk_blast.smk"
