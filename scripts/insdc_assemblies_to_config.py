#!/usr/bin/env python3

"""
Usage:
  insdc_assemblies_to_config.py <config_file>
    [--nodes /path/to/nodes.dmp] [--rank genus] [--root 2759] [--strategy WGS]
    [--out /path/to/output/directory]

Options:
  --nodes=<nodes>        Path to NCBI nodes dump file [default: nodes.dmp]
  --rank=<rank>          Similarity search database masking level
  --root=<root>          Root taxID [default: 2759] (default is all Eukaryota)
  --out=<out>            Path to output directory [default: .]
  --strategy=<strategy>  Sequencing strategy [default: WGS]
"""

import os
import re
import requests
import sys
import yaml
from docopt import docopt
from defusedxml import ElementTree as ET
from collections import defaultdict
from copy import deepcopy
import taxonomy

STEP = 50 # how many assembly records to request at a time

opts = docopt(__doc__)

NODES = opts['--nodes']
if not os.path.isfile(NODES):
    print("ERROR: File '%s' does not exist" % NODES, file=sys.stderr)
    quit(__doc__)

RANK = opts['--rank']
# awk '{print $5}' nodes.dmp | grep -v no | sort | uniq
valid_ranks = ( 'class','cohort','family','forma','genus','infraclass',
                'infraorder','kingdom','order','parvorder','phylum','species',
                'subclass','subfamily','subgenus','subkingdom','suborder',
                'subphylum','subspecies','subtribe','superclass','superfamily',
                'superkingdom','superorder','superphylum','tribe','varietas')
if RANK and RANK not in valid_ranks:
    print("ERROR: '%s' is not a valid value for --rank" % RANK, file=sys.stderr)
    quit(__doc__)

try:
    ROOT = int(opts['--root'])
except:
    print("ERROR: '%s' is not a valid value for --root" % opts['--root'], file=sys.stderr)
    quit(__doc__)

STRATEGY = opts['--strategy']

OUTDIR = "%s/%s" % (opts['--out'],ROOT)

try:
    os.makedirs(OUTDIR, exist_ok=True)
except:
    print("ERROR: Unable to create output directory '%s'" % OUTDIR, file=sys.stderr)
    quit(__doc__)

DEFAULT_META = {}

WITH_READS = "%s/sra" % OUTDIR
WITHOUT_READS = "%s/no_sra" % OUTDIR

os.makedirs("%s/" % WITH_READS, exist_ok=True)
os.makedirs("%s/" % WITHOUT_READS, exist_ok=True)

if os.path.isfile(sys.argv[1]):
    with open(sys.argv[1], 'r') as fh:
        try:
            DEFAULT_META = yaml.load(fh)
        except yaml.YAMLError as exc:
            print(exc)

"""
Generate config files for all INSDC assemblies from a given <ROOT> taxon for
which read data are available.
Group generated files in directories by <RANK> to simplify reuse of filtered
BLAST databases.
"""

def count_assemblies(root):
    """
    Query INSDC assemblies descended from <root> taxon.
    Return count as int.
    """
    warehouse = 'https://www.ebi.ac.uk/ena/data/warehouse'
    url = ("%s/search?query=\"tax_tree(%d)\"&result=assembly&resultcount"
        % (warehouse,root))
    response = requests.get(url)
    if response.ok:
        m = re.search('[\d,+]+',response.content.decode('utf-8'))
        return int(m.group(0).replace(',', ''))
    else:
        return 0

def deep_find_text(data,tags):
    """
    Find nested attributes in xml.
    Return attribute value.
    """
    for tag in tags:
        try:
            data = data.find(tag)
        except:
            return None
    return data.text

def list_assemblies(root,offset,count):
    """
    Query INSDC assemblies descended from <root> taxon.
    Return list of <count> entries from <offset> as tree.
    """
    warehouse = 'https://www.ebi.ac.uk/ena/data/warehouse'
    url = ("%s/search?query=\"tax_tree(%d)\"&result=assembly&display=xml&offset=%d&length=%d"
        % (warehouse,root,offset,count))
    response = requests.get(url)
    if response.ok:
        return response.content
    else:
        return None

def assembly_meta(asm,default_meta):
    """
    Return dict of metadata values for an assembly
    """
    meta = deepcopy(default_meta)
    if 'assembly' not in meta:
        meta['assembly'] = {}
    if 'taxon' not in meta:
        meta['taxon'] = {}
    genome_representation = asm.find('GENOME_REPRESENTATION').text
    if genome_representation == 'full':
        meta['assembly']['accession'] = asm.attrib['accession']
        meta['assembly']['bioproject'] = deep_find_text(asm,('STUDY_REF','IDENTIFIERS','PRIMARY_ID'))
        meta['assembly']['biosample'] = deep_find_text(asm,('SAMPLE_REF','IDENTIFIERS','PRIMARY_ID'))
        meta['taxon']['taxid'] = int(deep_find_text(asm,('TAXON','TAXON_ID')))
        meta['taxon']['name'] = deep_find_text(asm,('TAXON','SCIENTIFIC_NAME'))
        meta['assembly']['level'] = asm.find('ASSEMBLY_LEVEL').text
        meta['assembly']['alias'] = asm.attrib['alias']
        wgs_prefix = deep_find_text(asm,('WGS_SET','PREFIX'))
        wgs_version = deep_find_text(asm,('WGS_SET','VERSION'))
        if wgs_prefix and wgs_version:
            meta['assembly']['prefix'] = "%s%s" % (wgs_prefix,wgs_version.zfill(2))
        attributes = asm.find('ASSEMBLY_ATTRIBUTES')
        for attribute in attributes.findall('ASSEMBLY_ATTRIBUTE'):
            if attribute.find('TAG').text == 'total-length':
                meta['assembly']['span'] = int(attribute.find('VALUE').text)
            elif attribute.find('TAG').text == 'scaffold-count':
                meta['assembly']['scaffold-count'] = int(attribute.find('VALUE').text)
    return meta

def assembly_reads(biosample):
    """
    Query INSDC reads for a <biosample>.
    Return a dict of SRA accession, FASTQ ftp url, md5 and file size.
    """
    warehouse = 'https://www.ebi.ac.uk/ena/data/warehouse'
    url = ("%s/filereport?accession=%s&result=read_run&fields=run_accession,fastq_bytes,base_count,library_strategy,library_selection,library_layout,instrument_platform,fastq_ftp"
        % (warehouse,biosample))
    response = requests.get(url)
    sra = None
    if response.ok:
        lines = response.content.decode('utf-8').splitlines()
        if len(lines) > 1:
            sra = []
            header = lines[0].split('\t')
            for line in lines[1:]:
                fields = line.split('\t')
                values = {}
                reads = False
                strat = False
                for i in range(0,len(header)):
                    value = fields[i]
                    if header[i] == 'fastq_bytes':
                        value = fields[i].split(';')
                        if int(value[0] or 0) > 0:
                            reads = True
                            if len(value) >= 2:
                                values.update({'library_layout':'PAIRED'})
                            else:
                                values.update({'library_layout':'SINGLE'})
                    if header[i] == 'library_strategy':
                        if value == STRATEGY:
                            strat = True
                    values.update({header[i]:value})
                if reads and strat:
                    if 'base_count' not in values:
                        values['base_count'] = [0]
                    sra.append(values)
    return sra

def ncbi_assembly_reads(biosample):
    """
    Query NCBI INSDC reads for a <biosample>.
    Return a dict of SRA accession, etc.
    """
    if biosample != 'SAMN01914755':
        return None
    eutils = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
    esearch = ("%s/esearch.fcgi?db=sra&term=%s&usehistory=y"
        % (eutils,biosample))
    response = requests.get(esearch)
    sra = None
    if response.ok:
        search = ET.fromstring(response.content)
        count = int(search.find('Count').text)
        if count > 0:
            print(biosample)
            webenv=search.find('WebEnv').text
            key=search.find('QueryKey').text
            #print(webenv)
            esummary = ("%s/esummary.fcgi?db=sra&query_key=%s&WebEnv=%s"
                % (eutils,key,webenv))
            #print(esummary)
            sum_response = requests.get(esummary)
            if sum_response.ok:
                summaries = ET.fromstring(sum_response.content).findall('DocSum')
                for docsum in summaries:
                    text=docsum.find('Item').text
                    #print(text)
            #quit()
    return sra


graph = taxonomy.node_graph(NODES)
parents = taxonomy.parents_at_rank(graph,ROOT,RANK)
asm_count = count_assemblies(ROOT)

def base_count(x):
    if isinstance(x['base_count'], list):
        return int(x['base_count'][0] or 0)
    else:
        return 0

step = STEP
for offset in range(0,asm_count+1,step):
    count = step if offset + step < asm_count else asm_count - offset + 1
    print("%d: %d" % (offset,count))
    xml = list_assemblies(ROOT,offset,count)
    assemblies = ET.fromstring(xml)
    for assembly in assemblies:
        meta = {}
        meta = assembly_meta(assembly,DEFAULT_META)
        if 'prefix' in meta['assembly'] and 'biosample' in meta['assembly'] and meta['assembly']['biosample']:
            if RANK:
                if str(meta['taxon']['taxid']) in parents:
                    meta['taxon'][RANK] = int(parents[str(meta['taxon']['taxid'])])
                    meta['similarity']['defaults']['mask_ids'] = [meta['taxon'][RANK]]
            if 'reads' not in meta:
                meta['reads'] = {}
            sra = assembly_reads(meta['assembly']['biosample'])
            base_counts = {}
            fastq_ftp = {}
            if not sra:
                sra = ncbi_assembly_reads(meta['assembly']['biosample'])
            if sra:
                print(meta['assembly']['prefix'])
                sra.sort(key=lambda x: base_count(x), reverse=True)
                platforms = defaultdict(dict)
                for d in sra:
                    platforms[d['instrument_platform']].update({d['run_accession']:d})
                for platform,data in platforms.items():
                    meta['reads'][platform] = {}
                    strategies = defaultdict(list)
                    for key, value in data.items():
                        strategies[value['library_strategy']].append(key)
                    for strategy,accessions in strategies.items():
                        paired = []
                        single = []
                        for acc in accessions:
                            if data[acc]['library_layout'] == 'PAIRED':
                                paired.append(acc)
                            else:
                                single.append(acc)
                            base_counts[acc] = int(data[acc]['base_count'] or 0)
                            fastq_ftp[acc] = data[acc]['fastq_ftp']
                        if paired or single:
                            meta['reads'][platform][strategy] = {}
                            if paired:
                                meta['reads'][platform][strategy]['paired'] = paired
                            if single:
                                meta['reads'][platform][strategy]['single'] = single
                short_n = 3
                long_n = 10
                strategy = 'WGS'
                meta['reads']['paired'] = []
                meta['reads']['single'] = []
                for platform in ('ILLUMINA','LS454'):
                    if platform in meta['reads']:
                        if 'paired' in meta['reads'][platform][strategy]:
                            new_reads = meta['reads'][platform][strategy]['paired'][:short_n]
                            meta['reads']['paired'].extend([acc,platform,base_counts[acc],fastq_ftp[acc]] for acc in new_reads)
                            short_n -= len(meta['reads']['paired'])
                        if short_n > 0:
                            if 'single' in meta['reads'][platform][strategy]:
                                new_reads = meta['reads'][platform][strategy]['single'][:short_n]
                                meta['reads']['single'].extend([acc,platform,base_counts[acc],fastq_ftp[acc]] for acc in new_reads)
                                short_n -= len(meta['reads']['single'])
                for platform in ('PACBIO_SMRT','OXFORD_NANOPORE'):
                    if platform in meta['reads']:
                        if 'single' in meta['reads'][platform][strategy]:
                            new_reads = meta['reads'][platform][strategy]['single'][:long_n]
                            meta['reads']['single'] = [[acc,platform,base_counts[acc],fastq_ftp[acc]] for acc in new_reads]
                with open("%s/%s.yaml" % (WITH_READS,meta['assembly']['prefix']), 'w') as fh:
                    fh.write(yaml.dump(meta))
            else:
                with open("%s/%s.yaml" % (WITHOUT_READS,meta['assembly']['prefix']), 'w') as fh:
                    fh.write(yaml.dump(meta))
