import os

"""
https://github.com/blobtoolkit/insdc-pipeline
https://blobtoolkit.genomehubs.org/pipeline/

Pipeline to run BlobTools on public assemblies
----------------------------------------------

Requirements:
 - BlobTools2 (https://github.com/blobtoolkit/blobtools2)
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)

Basic usage:
  snakemake -p --use-conda \
    --directory ~/workdir \
    --configfile example.yaml \
    -j 8

© 2018-19 Richard Challis (University of Edinburgh), MIT License
© 2019-20 Richard Challis (Wellcome Sanger Institute), MIT License
"""

singularity: "docker://genomehubs/blobtoolkit:latest"

include: 'scripts/functions.py'

check_config()

multicore = int(os.getenv('MULTICORE', 16))
maxcore = int(os.getenv('MAXCORE', 32))

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
        "%s.fasta" % asm,
        expand("%s.{sra}.bam.stats" % asm, sra=list_sra_accessions(reads)),
        expand("%s%s/{sra}_cov.json" % (asm, rev), sra=list_sra_accessions(reads)),
        "%s%s/%s_phylum_positions.json" % (asm, rev, config['similarity']['taxrule']),
        expand("%s.busco.{lineage}.tsv" % asm, lineage=config['busco']['lineages']),
        expand("%s%s/{lineage}_busco.json" % (asm, rev), lineage=config['busco']['lineages']),
        "%s%s/identifiers.json" % (asm, rev)


# fetch database files
include: 'rules/fetch_ncbi_db.smk'
include: 'rules/fetch_taxdump.smk'
include: 'rules/fetch_uniprot.smk'
include: 'rules/extract_uniprot.smk'
include: 'rules/make_uniprot_db.smk'
include: 'rules/fetch_busco_lineage.smk'
# fetch assembly files
include: 'rules/fetch_assembly.smk'
include: 'rules/fetch_fastq.smk'
include: 'rules/subsample_fastq.smk'
# map reads
include: 'rules/bwa_index.smk'
include: 'rules/map_reads.smk'
include: 'rules/bamtools_stats.smk'
# make filtered databases
include: 'rules/split_fasta.smk'
include: 'rules/make_taxid_list.smk'
include: 'rules/make_masked_lists.smk'
include: 'rules/make_diamond_db.smk'
# run similarity searches
include: 'rules/run_windowmasker.smk'
include: 'rules/run_blastn.smk'
include: 'rules/extract_nohit_sequences.smk'
include: 'rules/run_diamond_blastx.smk'
include: 'rules/unchunk_blast.smk'
# run busco
include: 'rules/run_busco.smk'
# run blobtools
include: 'rules/generate_metadata.smk'
include: 'rules/blobtoolkit_create.smk'
include: 'rules/blobtoolkit_add_hits.smk'
include: 'rules/blobtoolkit_add_cov.smk'
include: 'rules/blobtoolkit_add_busco.smk'
