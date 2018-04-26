#!/usr/bin/env python3

# usage:
#     ./insdc_assemblies_to_config.py default.yaml

import os
import re
import requests
import sys
import yaml
from defusedxml import ElementTree as ET
from collections import defaultdict
from copy import deepcopy
import taxonomy

# TODO: set these with command line options
NODES = 'nodes.dmp'
RANK = 'genus' # used to determine database masking level
ROOT = 2759 # Eukaryota
ROOT = 7088 # Lepidoptera
NODES = '/Users/rchallis/tmp/btk_taxonomy/nodes.dmp'
ROOT = 6231 # Nematoda
ROOT = 6157 # Platyhelminthes
ROOT = 7214 # Drosophilidae
ROOT = 6447 # Mollusca
OUTDIR = "/Users/rchallis/projects/blobtoolkit/snakemake/insdc-config/%d" % ROOT
STEP = 50 # how many assembly records to request at a time
DEFAULT_META = {}

WITH_READS = "%s/sra" % OUTDIR
WITHOUT_READS = "%s/no_sra" % OUTDIR

os.makedirs(os.path.dirname("%s/" % WITH_READS), exist_ok=True)
os.makedirs(os.path.dirname("%s/" % WITHOUT_READS), exist_ok=True)

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
    meta['assembly'] = {}
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
    url = ("%s/filereport?accession=%s&result=read_run&fields=run_accession,fastq_bytes,library_strategy,library_selection,library_layout,instrument_platform"
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
                for i in range(0,len(header)):
                    value = fields[i]
                    if header[i] == 'fastq_bytes':
                        value = fields[i].split(';')
                    values.update({header[i]:value})
                sra.append(values)
    return sra

graph = taxonomy.node_graph(NODES)
parents = taxonomy.parents_at_rank(graph,ROOT,RANK)
asm_count = count_assemblies(ROOT)

step = STEP
for offset in range(1,asm_count+1,step):
    count = step if offset + step < asm_count else asm_count - offset + 1
    print("%d: %d" % (offset,count))
    xml = list_assemblies(ROOT,offset,count)
    assemblies = ET.fromstring(xml)
    for assembly in assemblies:
        meta = {}
        meta = assembly_meta(assembly,DEFAULT_META)
        if 'prefix' in meta['assembly'] and 'biosample' in meta['assembly'] and meta['assembly']['biosample']:
            print(meta['assembly']['prefix'])
            meta['taxon'][RANK] = int(parents[str(meta['taxon']['taxid'])])
            if 'similarity' in meta and 'defaults' in meta['similarity']:
                meta['similarity']['defaults']['mask_ids'] = [meta['taxon'][RANK]]
            meta['reads'] = {}
            sra = assembly_reads(meta['assembly']['biosample'])
            if not sra:
                sra = assembly_reads(meta['assembly']['bioproject'])
            if sra:
                sra.sort(key=lambda x: int(x['fastq_bytes'][0] or 0), reverse=True)
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
                            # if sra[acc]['library_layout'] == 'PAIRED':
                            if len(data[acc]['fastq_bytes']) >= 2:
                                paired.append(acc)
                            else:
                                single.append(acc)
                        if paired or single:
                            meta['reads'][platform][strategy] = {}
                            if paired:
                                meta['reads'][platform][strategy]['paired'] = paired
                            if single:
                                meta['reads'][platform][strategy]['single'] = single
                with open("%s/%s.yaml" % (WITH_READS,meta['assembly']['prefix']), 'w') as fh:
                    fh.write(yaml.dump(meta))
            else:
                with open("%s/%s.yaml" % (WITHOUT_READS,meta['assembly']['prefix']), 'w') as fh:
                    fh.write(yaml.dump(meta))
