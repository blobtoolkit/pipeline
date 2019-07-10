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
        "%s.fasta" % asm,
        expand("%s.{sra}.bam.stats" % asm,sra=list_sra_accessions(reads)),
        expand("%s/{sra}_cov.json" % asm,sra=list_sra_accessions(reads)),
        "%s/%s_phylum_positions.json" % (asm,config['similarity']['taxrule']),
        expand("%s_{lineage}.tsv" % asm,lineage=config['busco']['lineages']),
        expand("%s/{lineage}_busco.json" % asm,lineage=config['busco']['lineages']),
        "%s/identifiers.json" % asm


include: 'rules/fetch_database_files.smk'
include: 'rules/make_filtered_databases.smk'
include: 'rules/fetch_assembly_files.smk'
include: 'rules/run_similarity_searches.smk'
include: 'rules/map_reads.smk'
include: 'rules/run_blobtools.smk'
include: 'rules/run_busco.smk'
