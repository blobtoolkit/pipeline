#!/usr/bin/env python3

"""
Summarise numbers of assemblies available.

Usage:
  available_assemblies.py
    [--taxdump /path/to/ncbi_new_taxdump] [--root 2759]
    [--strategy WGS]

Options:
  --taxdump=<taxdump>    Path to NCBI taxdump directory [default: ../taxdump]
  --root=<root>          Root taxID [default: 2759] (default is all Eukaryota)
  --strategy=<strategy>  Sequencing strategy [default: WGS]
"""

import os
import re
import requests
import sys
import yaml
from docopt import docopt
from pathlib import Path
from defusedxml import ElementTree as ET
from collections import defaultdict
from copy import deepcopy
import taxonomy

STEP = 100  # how many assembly records to request at a time

opts = docopt(__doc__)

NODES = Path(opts['--taxdump']) / 'nodes.dmp'
if not os.path.isfile(NODES):
    print("ERROR: File '%s' does not exist" % NODES, file=sys.stderr)
    quit(__doc__)

NAMES = Path(opts['--taxdump']) / 'names.dmp'
if not os.path.isfile(NAMES):
    print("ERROR: File '%s' does not exist" % NAMES, file=sys.stderr)
    quit(__doc__)

LINEAGES = Path(opts['--taxdump']) / 'taxidlineage.dmp'
if not os.path.isfile(LINEAGES):
    print("ERROR: File '%s' does not exist" % LINEAGES, file=sys.stderr)
    quit(__doc__)

try:
    ROOT = int(opts['--root'])
except ValueError:
    print("ERROR: '%s' is not a valid value for --root" % opts['--root'], file=sys.stderr)
    quit(__doc__)

STRATEGY = opts['--strategy']


def grep_line(path, taxid):
    """Return the first line in the file at path beginning with taxid."""
    file = open(path, "r")
    for line in file:
        if re.match(str(taxid)+r'\s', line):
            return line
    return False


def split_line(line):
    """Split a lineage into a list of ancestral taxids."""
    taxids = [int(x) for x in line.split()[2:-1]]
    taxids.reverse()
    return taxids


def count_assemblies(root):
    """
    Query INSDC assemblies descended from <root> taxon.

    Return count as int.
    """
    warehouse = 'https://www.ebi.ac.uk/ena/data/warehouse'
    url = ("%s/search?query=\"tax_tree(%d)\"&result=assembly&resultcount"
           % (warehouse, root))
    response = requests.get(url)
    if response.ok:
        m = re.search(r'[\d,+]+', response.content.decode('utf-8'))
        return int(m.group(0).replace(',', ''))
    else:
        return 0


def deep_find_text(data, tags):
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


def list_assemblies(root, offset, count):
    """
    Query INSDC assemblies descended from <root> taxon.

    Return list of <count> entries from <offset> as tree.
    """
    warehouse = 'https://www.ebi.ac.uk/ena/data/warehouse'
    url = ("%s/search?query=\"tax_tree(%d)\"&result=assembly&display=xml&offset=%d&length=%d"
           % (warehouse, root, offset, count))
    try:
        response = requests.get(url)
        if response.ok:
            return response.content
        else:
            return None
    except:
        result = list_assemblies(root, offset, count)
        return result


def assembly_meta(asm):
    """Return dict of metadata values for an assembly."""
    meta = {}
    genome_representation = asm.find('GENOME_REPRESENTATION').text
    if genome_representation == 'full':
        meta['accession'] = asm.attrib['accession']
        meta['biosample'] = deep_find_text(asm, ('SAMPLE_REF', 'IDENTIFIERS', 'PRIMARY_ID'))
        meta['taxid'] = int(deep_find_text(asm, ('TAXON', 'TAXON_ID')))
        meta['alias'] = asm.attrib['alias']
        wgs_prefix = deep_find_text(asm, ('WGS_SET', 'PREFIX'))
        wgs_version = deep_find_text(asm, ('WGS_SET', 'VERSION'))
        if wgs_prefix and wgs_version:
            meta['prefix'] = "%s%s" % (wgs_prefix, wgs_version.zfill(2))
        elif ' ' not in meta['alias']:
            meta['prefix'] = meta['alias'].replace('.','_')
        else:
            meta['prefix'] = meta['accession'].replace('.','_')
    return meta


def has_reads(biosample):
    """
    Query INSDC reads for a <biosample>.

    Return True or False.
    """
    warehouse = 'https://www.ebi.ac.uk/ena/data/warehouse'
    url = ("%s/filereport?accession=%s&result=read_run"
           % (warehouse, biosample))
    try:
        response = requests.get(url)
        sra = False
        if response.ok:
            lines = response.content.decode('utf-8').splitlines()
            if len(lines) > 1:
                sra=True
        return sra
    except:
        return False


asm_count = count_assemblies(ROOT)
graph = taxonomy.node_graph(NODES)
species = taxonomy.parents_at_rank(graph, ROOT, 'species')
kingdoms = taxonomy.parents_at_rank(graph, ROOT, 'kingdom')


def hosted_assemblies(string='all'):
    """Get prefixes of hosted datasets from BTK API."""
    btk = "https://blobtoolkit.genomehubs.org/api/v1/search/%s" % string
    response = requests.get(btk)
    hosted = set()
    if response.ok:
        data = yaml.full_load(response.text)
        for asm in data:
            hosted.add(asm['prefix'])
    return hosted


search_term = 'all'
with open(NAMES, 'r') as fh:
    lines = fh.readlines()
    for l in lines:
        l = l[:-3]
        parts = re.split(r'\t\|\t', l)
        if parts[0] == str(ROOT) and parts[3] == 'scientific name':
            search_term = parts[1]
live = hosted_assemblies(search_term)

table = defaultdict(dict)

available = defaultdict(lambda: {'assemblies': set(), 'species': set(), 'with': set(), 'without': set()})
hosted = defaultdict(lambda: {'assemblies': set(), 'species': set(), 'with': set(), 'without': set()})
step = STEP
kingdom_ids = set()
offset = 0
processed = 0
processed_ids = set()
previous_processed = -1
while processed < asm_count:
    if processed == previous_processed:
        break
    previous_processed = processed
    count = step if offset + step < asm_count else asm_count - offset + 1
    print("%d: %d" % (offset, offset+count), file=sys.stderr)
    xml = list_assemblies(ROOT, offset, count)
    assemblies = ET.fromstring(xml)
    for assembly in assemblies:
        meta = {}
        meta = assembly_meta(assembly)
        if 'prefix' in meta and 'biosample' in meta and meta['biosample']:
            if str(meta['taxid']) in kingdoms:
                kingdom_id = kingdoms[str(meta['taxid'])]
                kingdom_ids.add(kingdom_id)
            else:
                kingdom_id = ROOT
                kingdom_ids.add(kingdom_id)
            if str(meta['taxid']) in species:
                species_id = int(species[str(meta['taxid'])])
            else:
                species_id = meta['taxid']
            if meta['prefix'] not in processed_ids:
                processed_ids.add(meta['prefix'])
                processed += 1
            print("INFO: %s has taxid %s (species: %s, kingdom: %s)" % (str(meta['prefix']), str(meta['taxid']), str(species_id), str(kingdom_id)), file=sys.stderr)
            available[kingdom_id]['species'].add(species_id)
            available[kingdom_id]['assemblies'].add(meta['prefix'])
            sra = has_reads(meta['biosample'])
            if sra:
                available[kingdom_id]['with'].add(meta['prefix'])
            else:
                available[kingdom_id]['without'].add(meta['prefix'])
            if meta['prefix'] in live:
                hosted[kingdom_id]['species'].add(species_id)
                hosted[kingdom_id]['assemblies'].add(meta['prefix'])
                if sra:
                    hosted[kingdom_id]['with'].add(meta['prefix'])
                else:
                    hosted[kingdom_id]['without'].add(meta['prefix'])
    offset += step
for kingdom_id in kingdom_ids:
    table[kingdom_id]['species'] = {'hosted': len(hosted[kingdom_id]['species']), 'available': len(available[kingdom_id]['species'])}
    table[kingdom_id]['assemblies'] = {}
    table[kingdom_id]['assemblies']['total'] = {'hosted': len(hosted[kingdom_id]['assemblies']), 'available': len(available[kingdom_id]['assemblies'])}
    table[kingdom_id]['assemblies']['with'] = {'hosted': len(hosted[kingdom_id]['with']), 'available': len(available[kingdom_id]['with'])}
    table[kingdom_id]['assemblies']['without'] = {'hosted': len(hosted[kingdom_id]['without']), 'available': len(available[kingdom_id]['without'])}
print(yaml.dump(dict(table)))
