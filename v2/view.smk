"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to validate a BlobDir and generate static views
--------------------------------------------------------

Requirements:
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p \
    --directory ~/workdir \
    --configfile example.yaml \
    -s view.smk
    -j 1

© 2021 Richard Challis (Wellcome Sanger Institute), MIT License
"""

include: "scripts/functions.py"


rule all:
    """
    Dummy rule to set target of pipeline
    """
    input:
        "%s/CHECKSUM" % blobdir_name(config)


include: "rules/validate_dataset.smk"
include: "rules/generate_images.smk"
include: "rules/generate_summary.smk"
include: "rules/checksum_files.smk"
