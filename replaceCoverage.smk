import os

"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to run replace read coverage in a BlobDir
--------------------------------------------------

Requirements:
 - BlobTools2 (https://github.com/blobtoolkit/blobtools2)
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p --use-conda \
    --directory ~/workdir \
    --configfile example.yaml \
    -s replaceHits.smk
    -j 8

© 2018-19 Richard Challis (University of Edinburgh), MIT License
© 2019-20 Richard Challis (Wellcome Sanger Institute), MIT License
"""

singularity: "docker://genomehubs/blobtoolkit:1.1"

include: 'scripts/functions.py'

multicore = int(os.getenv('MULTICORE', 16))
maxcore = int(os.getenv('MAXCORE', 32))

use_singularity = check_config()

similarity = apply_similarity_search_defaults()
reads = get_read_info(config)
keep = False
if 'keep_intermediates' in config:
    keep = bool(config['keep_intermediates'])
asm = config['assembly']['prefix']
rev = ''
if 'revision' in config:
    if config['revision'] > 0:
        rev = '.'+str(config['revision'])


rule all:
    """
    Dummy rule to set blobDB as target of pipeline
    """
    input:
        expand("%s.{sra}.bam.stats" % asm, sra=list_sra_accessions(reads)),
        expand("%s/{sra}_cov.json" % asm, sra=list_sra_accessions(reads)),
        "%s.meta.updated" % asm


# fetch database files
include: 'rules/fetch_taxdump.smk'
# fetch assembly files
include: 'rules/fetch_fastq.smk'
include: 'rules/subsample_fastq.smk'
# map reads
include: 'rules/bwa_index.smk'
include: 'rules/map_reads.smk'
include: 'rules/bamtools_stats.smk'
# run blobtools
include: 'rules/generate_metadata.smk'
include: 'rules/blobtoolkit_add_cov.smk'
include: 'rules/blobtoolkit_add_meta.smk'
