import os

"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to generate sequence stats
--------------------------------------------------

Requirements:
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p \
    --directory ~/workdir \
    --configfile example.yaml \
    -s stats.smk
    -j 8

Author:
  Richard Challis

Contact:
  blobtoolkit@genomehubs.org

License:
  © 2021 Genome Research Limited  % config["assembly"]["prefix"], MIT License
"""

include: "scripts/functions.py"

minimap_path = "../minimap"
windowmasker_path = "../windowmasker"

rule all:
    """
    Dummy rule to define output
    """
    input:
        "%s.stats.gc.bed"  % config["assembly"]["prefix"],
        "%s.stats.gc_windows.bed"  % config["assembly"]["prefix"],
        "%s.stats.n.bed"  % config["assembly"]["prefix"],
        "%s.stats.n_windows.bed"  % config["assembly"]["prefix"],
        "%s.stats.masked.bed"  % config["assembly"]["prefix"],
        "%s.stats.masked_windows.bed"  % config["assembly"]["prefix"],
        "%s.stats.length.bed" % config["assembly"]["prefix"],
        expand("%s.stats.{sra}_cov.bed" % config["assembly"]["prefix"], sra=reads_by_prefix(config).keys()),
        expand("%s.stats.{sra}_cov_windows.bed" % config["assembly"]["prefix"], sra=reads_by_prefix(config).keys()),
        

include: "rules/run_mosdepth.smk"
include: "rules/get_chunked_stats.smk"
