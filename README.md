# SnakeMake pipeline to run BlobTools on public assemblies

Make a copy of config.yaml in your working directory and specify assembly, reads
and databases to run the following steps:

- fetch and prepare taxon filtered ncbi and uniprot databases for blast and diamond searches
- fetch public assemblies and associated reads by accession
- map reads to assembly
- run similarity searches
- create Blob

This is a very resource intensive pipeline wit individual steps requiring up to 128GB RAM. It can be run using a command similar to:
```
snakemake -p --use-conda --directory working_directory/ --resources tmpdir=128 -j 64
```
or on a GridEngine cluster:
```
snakemake -p --use-conda --directory working_directory/ -j 64 --resources tmpdir=128 --cluster "qsub -v PATH=$PATH -pe smp {threads}"
```

## Requirements
- [Conda](https://conda.io/docs/commands/conda-install.html)
- [BlobTools](https://github.com/DRL/blobtools)
- [SnakeMake](http://snakemake.readthedocs.io/en/stable/)
