#!/usr/bin/env python3

"""
Generate config files for BlobToolKit pipeline.

Usage:
  generate_config.py <ACCESSION>
    [--rank genus] [--root 2759] [--coverage 20]
    [--out /path/to/output/directory] [--db /path/to/database/directory]

Options:
  --rank=<rank>          Similarity search database masking level [default: genus]
  --root=<root>          Root taxID [default: 2759] (default is all Eukaryota)
  --coverage=INT         Maximum coverage for read mapping [default: 20]
  --out=<out>            Path to output directory [default: .]
  --db=<db>              Path to database directory [default: .]
"""

import os
import re
import subprocess
import sys
from operator import itemgetter
from subprocess import PIPE, Popen

import ujson
from defusedxml import ElementTree as ET
from docopt import docopt
from tolkein import tofetch, tofile, tolog

LOGGER = tolog.logger(__name__)

GOAT_API = "https://goat.genomehubs.org/api/v0.0.1"
ENA_API = "https://www.ebi.ac.uk/ena/browser/api"
GCA_NAME = re.compile(r"GCA_.+")
BUSCO_URL = "https://busco-data.ezlab.org/v5/data/lineages"


def find_busco_lineages(ancestors):
    """Work out which BUSCO sets to run for a given taxid."""
    LOGGER.info("Identifying relevant BUSCO lineages")
    BUSCO_SETS = {
        "422676": "aconoidasida",
        "7898": "actinopterygii",
        "5338": "agaricales",
        "155619": "agaricomycetes",
        "33630": "alveolata",
        "5794": "apicomplexa",
        "6854": "arachnida",
        "6656": "arthropoda",
        "4890": "ascomycota",
        "8782": "aves",
        "5204": "basidiomycota",
        "68889": "boletales",
        "3699": "brassicales",
        "134362": "capnodiales",
        "33554": "carnivora",
        "91561": "cetartiodactyla",
        "34395": "chaetothyriales",
        "3041": "chlorophyta",
        "5796": "coccidia",
        "28738": "cyprinodontiformes",
        "7147": "diptera",
        "147541": "dothideomycetes",
        "3193": "embryophyta",
        "33392": "endopterygota",
        "314146": "euarchontoglires",
        "33682": "euglenozoa",
        "2759": "eukaryota",
        "5042": "eurotiales",
        "147545": "eurotiomycetes",
        "9347": "eutheria",
        "72025": "fabales",
        "4751": "fungi",
        "314147": "glires",
        "1028384": "glomerellales",
        "5178": "helotiales",
        "7524": "hemiptera",
        "7399": "hymenoptera",
        "5125": "hypocreales",
        "50557": "insecta",
        "314145": "laurasiatheria",
        "147548": "leotiomycetes",
        "7088": "lepidoptera",
        "4447": "liliopsida",
        "40674": "mammalia",
        "33208": "metazoa",
        "6029": "microsporidia",
        "6447": "mollusca",
        "4827": "mucorales",
        "1913637": "mucoromycota",
        "6231": "nematoda",
        "33183": "onygenales",
        "9126": "passeriformes",
        "5820": "plasmodium",
        "92860": "pleosporales",
        "38820": "poales",
        "5303": "polyporales",
        "9443": "primates",
        "4891": "saccharomycetes",
        "8457": "sauropsida",
        "4069": "solanales",
        "147550": "sordariomycetes",
        "33634": "stramenopiles",
        "32523": "tetrapoda",
        "155616": "tremellomycetes",
        "7742": "vertebrata",
        "33090": "viridiplantae",
        "71240": "eudicots",
    }
    # line = grep_line(lineage_file, taxid)
    lineages = []
    for obj in ancestors:
        if obj["taxon_id"] in BUSCO_SETS:
            lineages.append("%s_odb10" % BUSCO_SETS[obj["taxon_id"]])
    return lineages


def fetch_assembly_url(accession):
    """
    Fetch an assembly url using edirect.
    """
    esearch = Popen(
        "esearch -db assembly -query %s" % accession, stdout=PIPE, shell=True
    )
    esummary = Popen("esummary", stdin=esearch.stdout, stdout=PIPE, shell=True)
    xtract = Popen(
        "xtract -pattern DocumentSummary -element FtpPath_GenBank",
        stdin=esummary.stdout,
        stdout=PIPE,
        shell=True,
    )
    xtract_stdout = xtract.communicate()[0].decode("utf-8").strip()
    for url in xtract_stdout.split("\n"):
        basename = re.findall(GCA_NAME, url)[0]
        return "%s/%s_genomic.fna.gz" % (url, basename)
    return None


def fetch_assembly_fasta(url, filename):
    """Save assembly fasta file to local disk."""
    LOGGER.info("Fetching assembly FASTA to %s" % filename)
    tofetch.fetch_ftp(url, filename)


def fetch_assembly_meta_xml(accession):
    """
    Fetch assembly metadata xml from ENA.
    """
    url = "%s/xml/%s" % (ENA_API, accession)
    xml = tofetch.fetch_url(url)
    return xml


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


def parse_assembly_meta(accession):
    """Return dict of metadata values for an assembly."""
    LOGGER.info("Fetching assembly metadata")
    meta = {
        "assembly": {"accession": accession},
        "busco": {"lineages": []},
        "reads": {"paired": []},
        "settings": {
            "blast_chunk": 100000,
            "blast_max_chunks": 10,
            "blast_overlap": 500,
            "chunk": 1000000,
            "tmp": "/tmp",
        },
        "similarity": {
            "defaults": {
                "evalue": 1e-25,
                "mask_ids": [],
                "max_target_seqs": 10,
                "root": 1,
            },
            "taxrule": "bestdist",
            "databases": [],
        },
        "taxon": {},
    }
    xml = fetch_assembly_meta_xml(accession)
    root = ET.fromstring(xml)
    asm = root.find("ASSEMBLY")
    meta["assembly"]["bioproject"] = deep_find_text(
        asm, ("STUDY_REF", "IDENTIFIERS", "PRIMARY_ID")
    )
    meta["assembly"]["biosample"] = deep_find_text(
        asm, ("SAMPLE_REF", "IDENTIFIERS", "PRIMARY_ID")
    )
    meta["taxon"]["taxid"] = deep_find_text(asm, ("TAXON", "TAXON_ID"))
    meta["taxon"]["name"] = deep_find_text(asm, ("TAXON", "SCIENTIFIC_NAME"))
    meta["assembly"]["level"] = asm.find("ASSEMBLY_LEVEL").text
    meta["assembly"]["alias"] = asm.attrib["alias"]
    wgs_prefix = deep_find_text(asm, ("WGS_SET", "PREFIX"))
    wgs_version = deep_find_text(asm, ("WGS_SET", "VERSION"))
    if wgs_prefix and wgs_version:
        meta["assembly"]["prefix"] = "%s%s" % (wgs_prefix, wgs_version.zfill(2))
    elif " " not in meta["assembly"]["alias"]:
        meta["assembly"]["prefix"] = meta["assembly"]["alias"].replace(".", "_")
    else:
        meta["assembly"]["prefix"] = meta["assembly"]["accession"].replace(".", "_")
    attributes = asm.find("ASSEMBLY_ATTRIBUTES")
    for attribute in attributes.findall("ASSEMBLY_ATTRIBUTE"):
        if attribute.find("TAG").text == "total-length":
            meta["assembly"]["span"] = int(attribute.find("VALUE").text)
        elif attribute.find("TAG").text == "scaffold-count":
            meta["assembly"]["scaffold-count"] = int(attribute.find("VALUE").text)
    return meta


def fetch_busco_lineages(busco_sets, buscodir):
    """Fetch busco lineages."""
    if not busco_sets:
        return
    lineages_to_fetch = []
    for lineage in busco_sets:
        busco_lineage = "%s/%s" % (buscodir, lineage)
        if not os.path.isdir(busco_lineage):
            lineages_to_fetch.append(lineage)
    if not lineages_to_fetch:
        return
    lineage_urls = {}
    LOGGER.info("Fetching BUSCO lineage directory listing")
    listing = tofetch.fetch_url("%s/" % BUSCO_URL)
    for entry in listing.split("\n"):
        parts = re.split(r"[\"\s]+", entry)
        if len(parts) == 8:
            busco_set = re.sub(r"\..+$", "", parts[2])
            lineage_urls.update({busco_set: "%s/%s" % (BUSCO_URL, parts[2])})
    for lineage in lineages_to_fetch:
        LOGGER.info("Fetching BUSCO lineage %s" % lineage)
        tofetch.fetch_tar(lineage_urls[lineage], buscodir)


def fetch_goat_data(taxon_id):
    """Fetch taxon metadata from GoaT."""
    LOGGER.info("Fetching taxon metadata")
    url = "%s/record?recordId=taxon_id-%s&result=taxon" % (GOAT_API, taxon_id)
    result = tofetch.fetch_url(url)
    if result is None:
        LOGGER.error("Unable to fetch taxon metadata for '%s' from GoaT", taxon_id)
        sys.exit(1)
    data = ujson.loads(result)
    return data["records"][0]["record"]


def assembly_reads(biosample):
    """
    Query INSDC reads for a <biosample>.

    Return a dict of SRA accession, FASTQ ftp url, md5 and file size.
    """
    warehouse = "https://www.ebi.ac.uk/ena/data/warehouse"
    url = (
        "%s/filereport?accession=%s&result=read_run&fields=run_accession,fastq_bytes,base_count,library_strategy,library_selection,library_layout,instrument_platform,experiment_title,fastq_ftp"
        % (warehouse, biosample)
    )
    data = tofetch.fetch_url(url)
    sra = []
    header = None
    for line in data.split("\n"):
        if not line or line == "":
            continue
        if header is None:
            header = line.split("\t")
            continue
        fields = line.split("\t")
        values = {}
        x_ten = False
        for i in range(0, len(header)):
            value = fields[i]
            if header[i] == "experiment_title":
                if value == "HiSeq X Ten paired end sequencing":
                    x_ten = True
            values.update({header[i]: value})
        if x_ten:
            if "base_count" in values:
                values["base_count"] = int(values["base_count"])
            else:
                values["base_count"] = 0
            sra.append(values)
    if sra:
        return sorted(sra, key=itemgetter("base_count"), reverse=True)
    return None


def base_count(x):
    """Return number of bases or zero."""
    if isinstance(x["base_count"], list):
        return int(x["base_count"][0] or 0)
    else:
        return 0


def fetch_subsampled_reads(meta, ratio, readdir):
    """Fetch sra reads, subsampling if necessary."""
    for index, url in enumerate(meta["fastq_ftp"].split(";")):
        url = "ftp://%s" % url
        read_file = "%s/%s_%d.fastq.gz" % (readdir, meta["run_accession"], index + 1)
        if ratio > 1:
            LOGGER.info("Fetching read file %s", read_file)
            tofetch.fetch_ftp(url, read_file)
        else:
            LOGGER.info("Fetching read file %s", read_file)
            tmp_file = read_file.replace(".gz", ".tmp.gz")
            tofetch.fetch_ftp(url, tmp_file)
            LOGGER.info("Subsampling reads at ratio %.2f", ratio)
            cmd = "seqtk sample -s 100 %s %.2f" % (tmp_file, ratio)
            seqtk = Popen(cmd, stdout=PIPE, shell=True)
            gzip = Popen("gzip -c > %s" % read_file, stdin=seqtk.stdout, shell=True)
            gzip.wait()
            os.remove(tmp_file)


def add_taxon_to_meta(meta, taxon_meta):
    """Add taxon info to metadata."""
    LOGGER.info("Adding taxon metadata to assembly metadata")
    ranks = [
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "kingdom",
        "superkingdom",
    ]
    for obj in taxon_meta["lineage"]:
        if obj["taxon_rank"] in ranks:
            meta["taxon"].update({obj["taxon_rank"]: obj["scientific_name"]})
        if obj["taxon_rank"] == "genus":
            meta["similarity"]["defaults"].update({"mask_ids": [int(obj["taxon_id"])]})


def add_reads_to_meta(meta, sra, readdir):
    """Add read accessions to metadata."""
    LOGGER.info("Adding read accessions to assembly metadata")
    for index, run in enumerate(sra):
        info = [
            run["run_accession"],
            run["instrument_platform"],
            run["base_count"],
            "%s/%s_1.fastq.gz;%s/%s_2.fastq.gz"
            % (readdir, run["run_accession"], readdir, run["run_accession"]),
        ]
        meta["reads"]["paired"].append(info)
        if index == 2:
            return


# def current_versions(string="all"):
#     """Get curent versions of hosted datasets from BTK API."""
#     btk = "https://blobtoolkit.genomehubs.org/api/v1/search/%s" % string
#     response = requests.get(btk)
#     current = {}
#     if response.ok:
#         data = yaml.full_load(response.text)
#         for asm in data:
#             if "version" in asm:
#                 current.update({asm["prefix"]: asm["version"]})
#             else:
#                 current.update({asm["prefix"]: 1})
#     return current


# def create_outdir(span, version=1, _lineage="all"):
#     """Create output directory."""
#     span = tolkein.tobin.readable_bin(span)
#     name = "%s/v%s/%s" % (OUTDIR, str(version), span)
#     os.makedirs("%s/" % name, exist_ok=True)
#     return name


if __name__ == "__main__":
    opts = docopt(__doc__)
    accession = opts["<ACCESSION>"]
    outdir = opts["--out"]
    dbdir = opts["--db"]
    buscodir = "%s/busco" % dbdir
    uniprotdir = "%s/uniprot" % dbdir
    os.makedirs(buscodir, exist_ok=True)
    if not outdir.endswith(accession):
        outdir += "/%s" % accession
    os.makedirs(outdir, exist_ok=True)
    assembly_url = fetch_assembly_url(accession)
    os.makedirs("%s/assembly" % outdir, exist_ok=True)
    assembly_file = "%s/assembly/%s.fasta.gz" % (outdir, accession)
    fetch_assembly_fasta(assembly_url, assembly_file)
    meta = parse_assembly_meta(accession)
    meta["assembly"].update({"file": assembly_file})
    taxon_meta = fetch_goat_data(meta["taxon"]["taxid"])
    add_taxon_to_meta(meta, taxon_meta)
    #   import btk --fasta
    #   cumulative plot
    busco_sets = find_busco_lineages(taxon_meta["lineage"])
    if busco_sets:
        meta["busco"].update({"lineage_dir": buscodir, "lineages": busco_sets})
    fetch_busco_lineages(busco_sets, buscodir)
    #   run busco
    #   import btk --busco
    #   snail plot
    sra = assembly_reads(meta["assembly"]["biosample"])
    if sra:
        readdir = "%s/reads" % outdir
        add_reads_to_meta(meta, sra, readdir)
        os.makedirs(readdir, exist_ok=True)
        for run in sra:
            ratio = (
                meta["assembly"]["span"] * int(opts["--coverage"]) / run["base_count"]
            )
            # fetch_subsampled_reads(run, ratio, "%s/reads" % outdir)
    #   map reads
    #   import btk --cov
    #   greyscale blob
    #   generate filtered diamond database
    meta["similarity"]["databases"].append(
        {"local": uniprotdir, "name": "reference_proteomes",}
    )
    tofile.write_file("%s/config.yaml" % outdir, meta)
    #   run diamond blastx
    #   import btk --hits
    #   color blob


# search_term = "all"
# with open(NAMES, "r") as fh:
#     lines = fh.readlines()
#     for l in lines:
#         l = l[:-3]
#         parts = re.split(r"\t\|\t", l)
#         if parts[0] == str(ROOT) and parts[3] == "scientific name":
#             search_term = parts[1]
# versions = current_versions(search_term)

# step = STEP
# for offset in range(OFFSET, asm_count + 1, step):
#     count = step if offset + step < asm_count else asm_count - offset + 1
#     print("%d: %d" % (offset, count))
#     xml = list_assemblies(ROOT, offset, count)
#     assemblies = ET.fromstring(xml)
#     for assembly in assemblies:
#         meta = {}
#         meta = assembly_meta(assembly, DEFAULT_META)
#         if "span" not in meta["assembly"] or meta["assembly"]["span"] < MIN_SPAN:
#             continue
#         if (
#             "prefix" in meta["assembly"]
#             and "biosample" in meta["assembly"]
#             and meta["assembly"]["biosample"]
#         ):
#             if "similarity" not in meta:
#                 meta["similarity"] = {}
#             if "default" not in meta["similarity"]:
#                 meta["similarity"]["defaults"] = {}
#             if RANK:
#                 if str(meta["taxon"]["taxid"]) in parents:
#                     meta["taxon"][RANK] = int(parents[str(meta["taxon"]["taxid"])])
#                     meta["similarity"]["defaults"]["mask_ids"] = [meta["taxon"][RANK]]
#             if "reads" not in meta:
#                 meta["reads"] = {}
#             if "busco" not in meta:
#                 meta["busco"] = {}
#             meta["busco"]["lineages"] = find_busco_lineages(
#                 meta["taxon"]["taxid"], LINEAGES
#             )
#             sra = assembly_reads(meta["assembly"]["biosample"])
#             base_counts = {}
#             fastq_ftp = {}
#             # if not sra:
#             #     sra = ncbi_assembly_reads(meta['assembly']['biosample'])
#             version = 1
#             if meta["assembly"]["prefix"] in versions:
#                 version = versions[meta["assembly"]["prefix"]] + 1
#                 continue
#             meta["version"] = version
#             lineage = "all"
#             if meta["busco"]["lineages"]:
#                 lineage = meta["busco"]["lineages"][0]
#             print(meta["assembly"]["prefix"])
#             if sra:
#                 sra.sort(key=lambda x: base_count(x), reverse=True)
#                 platforms = defaultdict(dict)
#                 for d in sra:
#                     platforms[d["instrument_platform"]].update({d["run_accession"]: d})
#                 for platform, data in platforms.items():
#                     meta["reads"][platform] = {}
#                     strategies = defaultdict(list)
#                     for key, value in data.items():
#                         strategies[value["library_strategy"]].append(key)
#                     for strategy, accessions in strategies.items():
#                         paired = []
#                         single = []
#                         for acc in accessions:
#                             if data[acc]["library_layout"] == "PAIRED":
#                                 paired.append(acc)
#                             else:
#                                 single.append(acc)
#                             base_counts[acc] = int(data[acc]["base_count"] or 0)
#                             fastq_ftp[acc] = data[acc]["fastq_ftp"]
#                         if paired or single:
#                             meta["reads"][platform][strategy] = {}
#                             if paired:
#                                 meta["reads"][platform][strategy]["paired"] = paired
#                             if single:
#                                 meta["reads"][platform][strategy]["single"] = single
#                 short_n = 3
#                 long_n = 10
#                 strategy = "WGS"
#                 meta["reads"]["paired"] = []
#                 meta["reads"]["single"] = []
#                 for platform in ("ILLUMINA", "LS454"):
#                     if platform in meta["reads"]:
#                         if "paired" in meta["reads"][platform][strategy]:
#                             new_reads = meta["reads"][platform][strategy]["paired"][
#                                 :short_n
#                             ]
#                             meta["reads"]["paired"].extend(
#                                 [acc, platform, base_counts[acc], fastq_ftp[acc]]
#                                 for acc in new_reads
#                             )
#                             short_n -= len(meta["reads"]["paired"])
#                         if short_n > 0:
#                             if "single" in meta["reads"][platform][strategy]:
#                                 new_reads = meta["reads"][platform][strategy]["single"][
#                                     :short_n
#                                 ]
#                                 meta["reads"]["single"].extend(
#                                     [acc, platform, base_counts[acc], fastq_ftp[acc]]
#                                     for acc in new_reads
#                                 )
#                                 short_n -= len(meta["reads"]["single"])
#                 for platform in ("PACBIO_SMRT", "OXFORD_NANOPORE"):
#                     if platform in meta["reads"]:
#                         if "single" in meta["reads"][platform][strategy]:
#                             new_reads = meta["reads"][platform][strategy]["single"][
#                                 :long_n
#                             ]
#                             meta["reads"]["single"] = [
#                                 [acc, platform, base_counts[acc], fastq_ftp[acc]]
#                                 for acc in new_reads
#                             ]
#                 outdir = create_outdir(meta["assembly"]["span"], version, lineage)
#                 with open(
#                     "%s/%s.yaml" % (outdir, meta["assembly"]["prefix"]), "w"
#                 ) as fh:
#                     fh.write(yaml.dump(meta))
#             else:
#                 outdir = create_outdir(meta["assembly"]["span"], version, lineage)
#                 with open(
#                     "%s/%s.yaml" % (outdir, meta["assembly"]["prefix"]), "w"
#                 ) as fh:
#                     fh.write(yaml.dump(meta))
