import os

"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to run Blobtools
--------------------------------------------------

Requirements:
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p \
    --directory ~/workdir \
    --configfile example.yaml \
    -s blobtools.smk
    -j 8

Author:
  Richard Challis

Contact:
  blobtoolkit@genomehubs.org

License:
  © 2021 Genome Research Limited, MIT License
"""

include: "scripts/functions.py"

busco_path = "../busco"
minimap_path = "../minimap"
stats_path = "../stats"
diamond_path = "../diamond"
blastn_path = "../blastn"
diamond_blastp_path = "../diamond_blastp"

rule all:
    """
    Dummy rule to define output
    """
    input:
        "%s/meta.json" % blobdir_name(config),
        "%s/identifiers.json" % blobdir_name(config),
        "%s/%s_phylum.json" % (blobdir_name(config), similarity_setting(config, "diamond_blastx", "taxrule")),
        "%s/buscogenes_phylum.json" % blobdir_name(config),
        expand("%s/{sra}_cov.json" % blobdir_name(config), sra=reads_by_prefix(config).keys()),
        "%s/%s_busco.json" % (blobdir_name(config), config['busco']['lineages'][0])
        

include: "rules/run_blobtools_create.smk"
include: "rules/run_blobtools_add.smk"
include: "rules/add_summary_to_metadata.smk"
