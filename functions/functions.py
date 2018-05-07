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
            accessions += list(map(lambda sra: sra[0],reads['paired']))
        if 'single' in reads:
            accessions += list(map(lambda sra: sra[0],reads['single']))
    return accessions

def cov_files_by_platform(reads,assembly,platform):
    """
    Return a list of coverage files for a given sequencing platform.
    """
    accessions = []
    if reads is not None:
        if 'paired' in reads:
            accessions += [sra[0] for sra in reads['paired'] if sra[1] == platform]
        if 'single' in reads:
            accessions += [sra[0] for sra in reads['single'] if sra[1] == platform]
    return list(map(lambda sra: "%s.%s.bam.cov" % (assembly,sra),accessions))

def platform_cov_files(reads,assembly):
    platforms = set()
    if reads is not None:
        for strategy in ['paired','single']:
            for sra in reads[strategy]:
                if sra[1] not in platforms:
                    platforms.add(sra[1])
    return list(map(lambda platform: "%s.%s.sum.cov" % (assembly,platform),platforms))
