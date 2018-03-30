# apply defaults to similarity search databases
similarity = {}
for key,value in config['similarity']['defaults'].items():
    for db in config['similarity']['databases']:
        if key not in db:
            db[key] = value
        similarity.update({db['name']:db})

# TODO: add selected reads to reads->paired/reads->single
platforms = [k for k,v in config['reads'].items() if 'WGS' in v]
top_reads = {}
reads = {'paired':[],'single':[]}
for p in platforms:
    if p == 'CAPILLARY':
        continue
    top_reads.update({p:{}})
    ctr = 0
    if 'paired' in config['reads'][p]['WGS']:
        top_reads[p]['paired'] = []
        if len(config['reads'][p]['WGS']['paired']) > 1:
            top_reads[p]['paired'] = config['reads'][p]['WGS']['paired'][:2]
            ctr = 2
        else:
            top_reads[p]['paired'] = config['reads'][p]['WGS']['paired']
            ctr = 1
        reads['paired'] += top_reads[p]['paired']
    if ctr < 2:
        if 'single' in config['reads'][p]['WGS']:
            top_reads[p]['single'] = config['reads'][p]['WGS']['single'][0]
            ctr += 1
            if ctr < 2 and len(config['reads'][p]['WGS']['single']) > 1:
                top_reads[p]['single'].append(config['reads'][p]['WGS']['single'][1])
            reads['single'] += top_reads[p]['single']


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
