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

blobtools_path = "../blobtools"


rule all:
    """
    Dummy rule to set target of pipeline
    """
    input:
        "%s/CHECKSUM" % blobdir_name(config)


rule copy_blobdir:
    """Copy a blobdir to the working directory."""
    input:
        blobdir = "%s/{blobdir}" % blobtools_path
    output:
        copied = temp("{blobdir}.copied"),
        cov = expand("{{blobdir}}/{sra}_cov.json", sra=reads_by_prefix(config).keys()),
        tax = expand("{{blobdir}}/{taxrule}_phylum_positions.json", taxrule=config["similarity"]["taxrule"]),
        busco = expand("{{blobdir}}/{lineage}_busco.json", lineage=config["busco"]["lineages"]),
        ids = "{blobdir}/identifiers.json"
    params:
        blobdir = lambda wc: wc.blobdir
    threads: 1
    log:
        "logs/{blobdir}/copy_blobdir.log"
    benchmark:
        "logs/{blobdir}/copy_blobdir.benchmark.txt"
    shell:
        """(mkdir -p {params.blobdir} \
        && rsync -av {input.blobdir}/ {params.blobdir}/ \
        && touch {output.copied}) 2> {log}"""



include: "rules/validate_dataset.smk"
include: "rules/generate_images.smk"
include: "rules/generate_summary.smk"
include: "rules/checksum_files.smk"
