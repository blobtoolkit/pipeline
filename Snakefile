include: 'functions/functions.py'

similarity = apply_similarity_search_defaults()
reads = select_read_accessions()

print(reads)

rule all:
    """
    Dummy rule to set blobDB as target of pipeline
    """
    input:
        "%s.blobDB.json" % config['assembly']['prefix']

include: 'rules/fetch_database_files.smk'
include: 'rules/make_filtered_databases.smk'
include: 'rules/fetch_assembly_files.smk'
include: 'rules/run_similarity_searches.smk'
include: 'rules/map_short_reads.smk'
include: 'rules/run_blobtools.smk'
