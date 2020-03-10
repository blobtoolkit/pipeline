import re
import sys
import math

BWA_INDEX = ['amb', 'ann', 'bwt', 'pac', 'sa']


def check_config():
    """
    Check required fields are present in config.
    """
    sections = [{'name': 'assembly',
                 'keys': ['accession', 'prefix', 'alias', 'span'],
                 'defaults': {'accession': 'draft', 'alias': '==prefix'}},
                {'name': 'busco',
                 'keys': ['lineage_dir', 'lineages'],
                 'defaults': {'lineage_dir': 'busco_lineages', 'lineages': []}},
                {'name': 'reads',
                 'keys': ['paired', 'single'],
                 'defaults': {'paired': [], 'single': []}},
                {'name': 'settings',
                 'keys': ['blast_chunk', 'blast_max_chunks', 'blast_overlap',
                          'blobtools2_path', 'chunk', 'taxonomy', 'tmp'],
                 'defaults': {'blast_chunk': 100000, 'blast_max_chunks': 10,
                              'blast_overlap': 500, 'chunk': 1000000,
                              'tmp': '/tmp'}},
                {'name': 'similarity',
                 'keys': ['databases', 'defaults', 'taxrule'],
                 'defaults': {'defaults': {'evalue': 1e-25, 'max_target_seqs': 10, 'root': 1},
                              'taxrule': 'bestsumorder'}},
                {'name': 'taxon',
                 'keys': ['name', 'taxid'],
                 'defaults': {}}]
    optional = ['busco', 'reads']
    for section in sections:
        if section['name'] not in config:
            if section['name'] in optional:
                print("INFO: optional section '%s' is not present in config file" % section['name'])
                config[section['name']] = {}
            else:
                quit("ERROR: config file must contain a '%s' section with keys '%s'" % (section['name'],
                                                                                        ', '.join(section['keys'])))
        for key in section['keys']:
            if key not in config[section['name']]:
                if key in section['defaults']:
                    value = section['defaults'][key]
                    if isinstance(value, str) and value.startswith('=='):
                        value = config[section['name']][value.replace('==', '')]
                    print("INFO: using default value for '%s.%s'" % (section['name'], key))
                    print(value)
                    config[section['name']][key] = value
                else:
                    quit("ERROR: config file section '%s' must contain '%s'" % (section['name'], key))
    # fill in additional database info
    for db in config['similarity']['databases']:
        if 'name' not in db or 'local' not in db:
            quit("ERROR: 'name' and 'local' must be specified for all databases")
        if db['name'] == 'nt' or db['name'] == 'nt_v5':
            db.update({'source': 'ncbi',
                       'tool': 'blast',
                       'type': 'nucl'})
        elif db['name'] == 'reference_proteomes':
            db.update({'source': 'uniprot',
                       'tool': 'diamond',
                       'type': 'prot'})
        else:
            print("INFO: only 'nt' and 'reference_proteomes' databases are supported, ignoring '%s'" % db['name'])
    if '--use-singularity' in sys.argv:
        return True
    return False


def apply_similarity_search_defaults():
    """
    Apply defaults to similarity search databases.
    """
    similarity = {}
    for key, value in config['similarity']['defaults'].items():
        for db in config['similarity']['databases']:
            if key not in db:
                db[key] = value
            similarity.update({db['name']: db})
            if db['name'].startswith('nt'):
                similarity.update({'blastdb': {'local': db['local']}})
    return similarity


def get_read_info(config):
    """
    Create dict of sequencing strategies, platforms and base count for reads.
    """
    reads = {}
    min = 0
    max = math.inf
    platforms = ('ILLUMINA', 'OXFORD_NANOPORE', 'PACBIO_SMRT', 'LS454')
    strategies = ('paired', 'single')
    if 'reads' not in config:
        return reads
    if 'coverage' in config['reads']:
        if 'min' in config['reads']['coverage']:
            min = config['reads']['coverage']['min']
        if 'max' in config['reads']['coverage']:
            max = config['reads']['coverage']['max']
    for strategy in strategies:
        for row in config['reads'][strategy]:
            accession = row[0]
            platform = row[1]
            if platform not in platforms:
                print("WARNING: platform %s is not recognised, must be one of %s" % (platform, platforms),
                      file=sys.stderr)
            try:
                bases = row[2]
                coverage = bases / config['assembly']['span']
            except:
                coverage = 10
            if strategy == 'paired':
                try:
                    url = re.split(',|;', row[3])
                    if len(reads) > 2:
                        reads = reads[-2:]
                except:
                    url = ["%s_1.fastq.gz" % accession, "%s_2.fastq.gz" % accession]
            else:
                try:
                    url = [row[3]]
                except:
                    url = ["%s.fastq.gz" % accession]
            if coverage >= min:
                reads[accession] = {'platform': platform, 'coverage': coverage, 'strategy': strategy, 'url': url}
                if coverage > max:
                    reads[accession]['subsample'] = max / coverage
                    print("WARNING: read file %s will be subsampled due to high coverage (%.2f > %.2f)" % (accession, coverage, max),
                          file=sys.stderr)
            else:
                print("WARNING: skipping read file %s due to low coverage (%.2f < %.2f)" % (accession, coverage, min),
                      file=sys.stderr)
    return reads


def ncbi_idmap(name):
    """
    Make a list of remote "accession2taxid" files to download.
    """
    url = 'ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid'
    db = similarity[name]
    return ' '.join(list(map(lambda x: "%s/%s.accession2taxid.gz" % (url, x), db['idmap'])))


def list_similarity_results(config):
    """
    Generate a list of output filenames for sequence similarity searches
    based on list of databases in "config['similarity']".
    """
    path = []
    for db in config['similarity']['databases']:
        # suffix = 'out' if db['tool'] == 'blast' else 'taxified.out'
        suffix = 'out'
        program = 'blastn' if db['type'] == 'nucl' else 'blastx' if db['tool'] == 'blast' else 'diamond'
        masked = ''
        if 'mask_ids' in db and isinstance(db['mask_ids'], (list, )):
            masked = "minus.%s" % '.'.join(str(mask) for mask in db['mask_ids'])
        else:
            masked = 'full'
        path.append("%s.%s.%s.root.%s.%s.%s" % (config['assembly']['prefix'],
                                                program,
                                                db['name'],
                                                db['root'],
                                                masked,
                                                suffix))
    return path


def blast_db_name(config):
    """
    Test whether _v5 should be appended to nt database name.
    """
    for db in config['similarity']['databases']:
        if db['name'].startswith('nt'):
            return db['name']
    return 'nt'


def blast_query_file(name, assembly):
    """
    Generate filename for filtered query file for similarity searches.
    """
    file = assembly
    if 'exclude_hits' in similarity[name]:
        for db in similarity[name]['exclude_hits']:
            file = db + '_filtered.' + file
    return file


def list_sra_accessions(reads):
    """
    Return a list SRA accessions.
    """
    accessions = []
    if reads is not None:
        accessions = reads.keys()
    return accessions


def generate_mapping_command(accession, reads):
    """
    Generate a read mapping command appropriate to the
    sequencing strategy and library type.
    """
    cmd = 'bwa mem'
    if reads[accession]['platform'] == 'ILLUMINA':
        cmd = 'minimap2 -ax sr'
    elif reads[accession]['platform'] == 'PACBIO_SMRT':
        cmd = 'minimap2 -ax map-pb'
    elif reads[accession]['platform'] == 'OXFORD_NANOPORE':
        cmd = 'minimap2 -ax map-ont'
    return cmd


def list_read_files(accession, reads, subsample):
    """
    List read files.
    """
    files = []
    for fq_url in reads[accession]['url']:
        file = re.sub(r'.+\/', '', fq_url)
        if subsample and 'subsample' in reads[accession]:
            file = file.replace('fastq', 'subsampled.fastq')
        files.append(file)
    return files


def generate_subsample_command(accession, reads):
    """
    Generate a read mapping command appropriate to the
    sequencing strategy and library type.
    """
    cmd = 'cp'
    arrow = ''
    seed = 100
    if 'coverage' in reads and 'seed' in reads['coverage']:
        seed = reads['coverage']['seed']
    if 'subsample' in reads[accession]:
        cmd = "seqtk sample -s%s" % seed
        arrow = "%.2f | pigz -c > " % reads[accession]['subsample']
    return [cmd, arrow]


def prepare_ebi_sra_url(acc, file):
    if len(reads[acc]) == 1:
        return reads[acc][0]
    urls = []
    for url in reads[acc]['url']:
        urls += url.split(';')
    for url in urls:
        if file in url:
            if url.startswith('ftp') and 'ftp://' not in url:
                url = 'ftp://'+url
            return url
    return ''


def prepare_ncbi_assembly_url(accession, name):
    base = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all'
    acc = accession.replace('_', '').split('.', 1)[0]
    path = '/'.join(acc[i:i+3] for i in range(0, len(acc), 3))
    asm = "%s_%s" % (accession, name.replace(' ', '_'))
    asm = asm.replace('__', '_').replace(',', '')
    url = "%s/%s/%s/%s_genomic.fna.gz" % (base, path, asm, asm)
    return url


def cov_files_by_platform(reads, assembly, platform):
    """
    Return a list of coverage files for a given sequencing platform.
    """
    accessions = []
    if reads is not None:
        accessions += [accession for accession in reads if reads[accession]['platform'] == platform]
    return list(map(lambda sra: "%s.%s.bam.cov" % (assembly, sra), accessions))


def platform_cov_files(reads, assembly):
    platforms = set()
    if reads is not None:
        for accession in reads:
            if reads[accession]['platform'] not in platforms:
                platforms.add(reads[accession]['platform'])
    return list(map(lambda platform: "%s.%s.sum.cov" % (assembly, platform), platforms))


def get_threads(rule, default, ratio=1):
    if rule in cluster_config and 'threads' in cluster_config[rule]:
        return int(cluster_config[rule]['threads'] * ratio)
    elif '__default__' in cluster_config and 'threads' in cluster_config['__default__']:
        return int(cluster_config['__default__']['threads'] * ratio)
    return default


def hit_fields(asm, rev, taxrule):
    if taxrule.startswith('best'):
        return "%s%s/%s_phylum_positions.json" % (asm, rev, taxrule)
    if taxrule.startswith('each'):
        return [
            "%s%s/%s_nt_phylum_positions.json" % (asm, rev, taxrule.replace('each', 'best')),
            "%s%s/%s_aa_phylum_positions.json" % (asm, rev, taxrule.replace('each', 'best'))
        ]


def ncbi_dir(wc):
    if use_singularity:
        dir = '/blobtoolkit/databases/ncbi_db'
    else:
        dir = similarity['blastdb']['local']
    return dir


def uniprot_dir(wc):
    if use_singularity:
        dir = '/blobtoolkit/databases/uniprot_db'
    else:
        dir = similarity['reference_proteomes']['local']
    return dir


def busco_dir(wc):
    if use_singularity:
        dir = '/blobtoolkit/databases/busco'
    else:
        dir = config['busco']['lineage_dir']
    return dir


def taxdump_dir(wc):
    if use_singularity:
        dir = '/blobtoolkit/databases/ncbi_taxdump'
    else:
        dir = config['settings']['taxonomy']
    return dir


def git_dir(wc):
    if use_singularity:
        dir = '/blobtoolkit/insdc-pipeline/.git'
    else:
        dir = os.path.dirname(os.path.abspath(workflow.snakefile))+'/.git'
    return dir
