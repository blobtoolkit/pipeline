BWA_INDEX = ['amb','ann','bwt','pac','sa']

def apply_similarity_search_defaults():
    """
    Apply defaults to similarity search databases.
    """
    similarity = {}
    if 'defaults' in config['similarity']:
        for key,value in config['similarity']['defaults'].items():
            for db in config['similarity']['databases']:
                if key not in db:
                    db[key] = value
                similarity.update({db['name']:db})
    else:
        for db in config['similarity']['databases']:
            similarity.update({db['name']:db})
    return similarity

def select_read_accessions(strategy):
    """
    Select read accessions to map to an assembly (up to 2 read accessions per
    platform) unless paired/single read accessions are already specified.
    Return a dict containing lists of read accessions to use.
    """
    if 'paired' in config['reads'] or 'single' in config['reads']:
        return config['reads']
    reads = {'paired':[],'single':[]}
    platforms = [k for k,v in config['reads'].items() if strategy in config['reads'][k]]
    if not platforms:
        reads = None
    top_reads = {}
    for p in platforms:
        top_reads.update({p:{}})
        ctr = 0
        if 'paired' in config['reads'][p][strategy]:
            top_reads[p]['paired'] = []
            if len(config['reads'][p][strategy]['paired']) > 1:
                top_reads[p]['paired'] = config['reads'][p][strategy]['paired'][:2]
                ctr = 2
            else:
                top_reads[p]['paired'] = config['reads'][p][strategy]['paired']
                ctr = 1
            reads['paired'] += top_reads[p]['paired']
        if ctr < 2:
            if 'single' in config['reads'][p][strategy]:
                top_reads[p]['single'] = [config['reads'][p][strategy]['single'][0]]
                ctr += 1
                if ctr < 2 and len(config['reads'][p][strategy]['single']) > 1:
                    top_reads[p]['single'].append(config['reads'][p][strategy]['single'][1])
                reads['single'] += top_reads[p]['single']
    return reads

def ncbi_idmap(name):
    """
    Make a list of remote "accession2taxid" files to download
    """
    url = 'ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid'
    db = similarity[name]
    return ' '.join(list(map(lambda x: "%s/%s.accession2taxid.gz" % (url,x),db['idmap'])))

def list_similarity_results(config):
    """
    Generate a list of output filenames for sequence similarity searches
    based on list of databases in "config['similarity']".
    """
    path = []
    for db in config['similarity']['databases']:
        suffix = 'out' if db['tool'] == 'blast' else 'taxified.out'
        program = 'blastn' if db['type'] == 'nucl' else 'blastx' if db['tool'] == 'blast' else 'diamond'
        masked = ''
        if 'mask_ids' in db and isinstance(db['mask_ids'],(list,)):
            masked = "minus.%s" % '.'.join(str(mask) for mask in db['mask_ids'])
        else:
            masked = 'full'
        path.append("%s.%s.%s.root.%s.%s.%s" % (config['assembly']['prefix'],program,db['name'],db['root'],masked,suffix))
    return path

def list_sra_accessions(reads):
    """
    Return a list SRA accessions.
    """
    accessions = []
    if reads is not None:
        if 'paired' in reads:
            accessions += list(map(lambda sra: sra,reads['paired']))
        if 'single' in reads:
            accessions += list(map(lambda sra: sra,reads['single']))
    return accessions
