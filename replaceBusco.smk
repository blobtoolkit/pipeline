"""
https://github.com/blobtoolkit/insdc-pipeline

Pipeline to run BlobTools on public assemblies
----------------------------------------------

Requirements:
 - BLAST+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
 - BlobTools2 (https://github.com/blobtoolkit/blobtools2)
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p --use-conda
    --directory path/to/workdir/
    --configfile path/to/config.yaml
    -j 8

© 2018-19 Richard Challis (University of Edinburgh), MIT License
"""


include: 'scripts/functions.py'

similarity = apply_similarity_search_defaults()
reads = get_read_info(config)
keep = False
if 'keep_intermediates' in config:
    keep = bool(config['keep_intermediates'])
asm = config['assembly']['prefix']

rule all:
    """
    Dummy rule to set blobDB as target of pipeline
    """
    input:
        'busco.replaced',
        expand("%s/{lineage}_busco.json" % asm,lineage=config['busco']['lineages'])


include: 'rules/blobtools_replace.smk'
