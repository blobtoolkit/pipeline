configfile: 'config.yaml'

ENV = config['settings']['environment']

# apply defaults to similarity search databases
similarity = {}
for key,value in config['similarity']['defaults'].items():
    for db in config['similarity']['databases']:
        if key not in db:
            db[key] = value
        similarity.update({db['name']:db})

rule all:
    """
    Dummy rule to set blobDB as target of pipeline
    """
    input:
        "%s.blobDB.json" % config['assembly']['name']

include: 'rules/fetch_database_files.smk'
include: 'rules/make_filtered_databases.smk'
include: 'rules/fetch_assembly_files.smk'
include: 'rules/run_similarity_searches.smk'
include: 'rules/map_short_reads.smk'
include: 'rules/run_blobtools.smk'
