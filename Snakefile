"""
https://github.com/blobtoolkit/insdc-pipeline

Pipeline to run BlobTools on public assemblies
----------------------------------------------

Requirements:
 - BlobTools (https://github.com/DRL/blobtools)
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - enaBrowserTools (https://github.com/enasequence/enaBrowserTools)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p --use-conda
    --directory path/to/workdir/
    --configfile path/to/config.yaml
    -j 8

© 2018 Richard Challis (University of Edinburgh), MIT License
"""


include: 'functions/functions.py'

similarity = apply_similarity_search_defaults()
reads = config['reads']

asm = config['assembly']['prefix']

rule all:
    """
    Dummy rule to set blobDB as target of pipeline
    """
    input:
        "%s.fasta" % asm,
        "%s.blobDB.json" % asm,
        expand("%s.{sra}.bam.stats" % asm,sra=list_sra_accessions(reads))


include: 'rules/fetch_database_files.smk'
include: 'rules/make_filtered_databases.smk'
include: 'rules/fetch_assembly_files.smk'
include: 'rules/run_similarity_searches.smk'
include: 'rules/map_short_reads.smk'
include: 'rules/run_blobtools.smk'
