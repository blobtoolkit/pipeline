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

© 2021 Richard Challis (Wellcome Sanger Institute), MIT License
"""

include: "scripts/functions.py"

rule all:
    """
    Dummy rule to define output
    """
    input:
        "%s/meta.json" % blobdir_name(config),
        "%s/identifiers.json" % blobdir_name(config),
        "%s/%s_phylum.json" % (blobdir_name(config), config["similarity"]["taxrule"]),
        expand("%s/{sra}_cov.json" % blobdir_name(config), sra=reads_by_prefix(config).keys()),
        expand("%s/{lineage}_busco.json" % blobdir_name(config), lineage=config['busco']['lineages'])
        

include: "rules/run_blobtools_create.smk"
include: "rules/unzip_assembly_fasta.smk"
include: "rules/run_bamtools_stats.smk"
include: "rules/add_summary_to_metadata.smk"
# include: "rules/make_diamond_db.smk"
# include: "rules/chunk_fasta.smk"
# include: "rules/run_diamond_blastx.smk"
# include: "rules/unchunk_blast.smk"
