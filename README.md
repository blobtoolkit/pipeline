# insdc-pipeline

_A [SnakeMake](http://snakemake.readthedocs.io/en/stable/) pipeline to run [BlobTools](https://github.com/DRL/blobtools) on public genome assemblies_

Make a copy of config.yaml in your working directory and specify assembly, reads
and databases to run the following steps:

- fetch and prepare taxon filtered [NCBI](https://www.ncbi.nlm.nih.gov) and [UniProt](http://www.uniprot.org) databases for [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and [Diamond](https://github.com/bbuchfink/diamond) sequence similarity searches
- fetch public assemblies and associated reads by accession using [enaBrowserTools](https://github.com/enasequence/enaBrowserTools)
- map reads to the assembly using [bwa](http://bio-bwa.sourceforge.net)
- run similarity searches
- combine input files to create a [BlobTools](https://github.com/DRL/blobtools) BlobDB

## Requirements
The following programs need to be installed and available in your $PATH:
- [BlobTools](https://github.com/DRL/blobtools)
- [Conda](https://conda.io/docs/commands/conda-install.html)
- [enaBrowserTools](https://github.com/enasequence/enaBrowserTools) (python3)
- [SnakeMake](http://snakemake.readthedocs.io/en/stable/)

All other required programs and libraries can be installed automatically via [Conda](https://conda.io/docs/commands/conda-install.html) when the pipeline is run.

## Usage

Make a copy of `config.yaml` in your working directory and edit the relevant details to fetch fasta and fastq files for your chosen assembly. Config files for specific taxa can be generated using the script `insdc_assemblies_to_config.py`.

To run all steps on a single machine (using up to 8 cores):

```
snakemake -p --use-conda \
    --directory path/to/workdir \
    --configfile path/to/config.yaml \
    -j 8
```

To submit individual steps as jobs to a gridengine cluster (using settings in cluster.yaml):

```
BTKDIR=/path/to/blobtoolkit/insdc-pipeline
WORKDIR=/path/to/workdir
ASSEMBLY=AIXA01
> $BTKDIR/$SET/$ASSEMBLY.qsub.log &&
snakemake -p \
    --use-conda --conda-prefix /path/to/.conda \
    --cluster-config cluster.yaml \
    --cluster "qsub \
        -v PATH=/path/to/blobtools:/path/to/enaBrowserTools/python3:$PATH \
        -pe {cluster.pe} {cluster.threads} \
        -j y -o $WORKDIR/$ASSEMBLY.qsub.log" \
    --jobname "{rulename}.{jobid}.sh" \
    --directory $WORKDIR \
    --configfile $WORKDIR/$ASSEMBLY.yaml \
    --stats $ASSEMBLY.snakemake.stats \
    -j 64
```

Note that individual steps may require up to 128 GB RAM.
