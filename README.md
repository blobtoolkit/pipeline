# insdc-pipeline

_A [Snakemake](http://snakemake.readthedocs.io/en/stable/) pipeline to run [BlobTools](https://github.com/DRL/blobtools) on public genome assemblies_

While designed for use on public genome assemblies and read files, this pipeline is also suitable for use with local assembly and read files. If the assemby and read files exist locally when the pipeline is run, the remote fetches will be skipped.

## Overview

The final output of this pipeline is a collection of JSON files ready to be viewed in the [BlobToolKit Viewer](http://blobtoolkit.genomehubs.org/view/) (Figure 1).

![Figure 1](http://blobtoolkit.genomehubs.org/wp-content/uploads/2018/08/figure1.png)

The pipeline (Figure 2) is implemented using Snakemake and will fetch all required database and assembly files then run BLAST/Diamond similarity searches and bwa/minimap2 read mapping and for processing with BlobTools.

![Figure 2](http://blobtoolkit.genomehubs.org/wp-content/uploads/2018/08/figure2.png)

## Installation

Most dependencies are handled by Snakemake. To get started you will need to install Conda, BLAST+ and BlobTools, then use Conda to install Snakemake.

### Conda

Follow instructions at [conda.io](https://conda.io/docs/user-guide/install/index.html)

### BLAST+

Download relevant binary from [ftp.ncbi.nlm.nih.gov](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.0alpha/)

### Blobtools

This pipeline requires an experimental view option currently only up to date in this [Blobtools fork](https://github.com/rjchallis/blobtools).

```
$ curl https://raw.githubusercontent.com/blobtoolkit/insdc-pipeline/master/envs/blobtools.yaml > blobtools.yaml
$ conda env create -n blobtools_env -f blobtools.yaml
$ conda activate blobtools_env
$ git clone https://github.com/rjchallis/blobtools
$ cd blobtools
$ ./install
$ conda deactivate
```

### Snakemake

```
$ conda create -n snake_env -c bioconda snakemake
```

## Usage

Make a copy of `config.yaml` matching your assembly name in your working directory and edit the relevant details to fetch fasta and fastq files for your chosen assembly. Change all file/directory paths to match your system and the locations of the BLAST+ and BlobTools executables. Config files for public assemblies can also be generated using the script `insdc_assemblies_to_config.py`.

### BDQP01.yaml

```
assembly:
  accession: GCA_002217835.1
  alias: Dobs_1.0
  bioproject: PRJDB4578
  biosample: SAMD00047210
  level: scaffold
  prefix: BDQP01
  scaffold-count: 1935
  span: 181868570
reads:
  paired:
  - [DRR055278, ILLUMINA]
  - [DRR055279, ILLUMINA]
  - [DRR055280, ILLUMINA]
settings:
  chunk: 1000000
  blast_chunk: 100000
  blast_overlap: 500
  blobtools_path: /software/blobtools
  blast_path: /software/blast/ncbi-blast-2.8.0+/bin
  taxonomy: /software/databases/ncbi_2018_05/taxonomy
  tmp: /tmp
similarity:
  databases:
  - idmap: [nucl_est, nucl_gb, nucl_gss, nucl_wgs, pdb]
    local: /software/databases/ncbi_2018_05
    name: nt
    source: ncbi
    tool: blast
    type: nucl
  - local: /software/databases/uniprot_2018_04
    max_target_seqs: 1
    name: reference_proteomes
    source: uniprot
    tool: diamond
    type: prot
  defaults:
    evalue: 1e-25
    mask_ids: [7215]
    max_target_seqs: 10
    root: 1
  taxrule: bestsumorder
taxon:
  name: Drosophila obscura
  taxid: 7282
  family: Drosophilidae
  phylum: Arthropoda
```


```
git clone https://github.com/blobtoolkit/insdc-pipeline
conda activate snake_env
BTKDIR=`pwd`/insdc-pipeline
cd $BTKDIR
CORES=32
WORKDIR=/path/to/working/directory
ASSEMBLY=BDQP01;
snakemake -p --use-conda --conda-prefix ~/.conda --directory $WORKDIR/ --configfile $WORKDIR/$ASSEMBLY.yaml --stats $ASSEMBLY.snakemake.stats -j $CORES --restart-times 2
```


Note that individual steps may require up to 128 GB RAM.
