# SnakeMake pipeline to run BlobTools on public assemblies

Make a copy of config.yaml in your working directory and specify assembly, reads
and databases to run the following steps:

- fetch and prepare taxon filtered ncbi and uniprot databases for blast and diamond searches
- fetch public assemblies and associated reads by accession
- map reads to assembly
- run similarity searches
- create Blob

run using a command similar to:
```
snakemake -p --use-conda --directory working_directory/ --resources tmpdir=128 -j 64
```

## Before you run the pipeline
Please note that some steps have yet to be optimised so this is likely to take several days to run and requires a machine with at least 128GB RAM.

## Requirements
- [Conda](https://conda.io/docs/commands/conda-install.html)
- [BlobTools](https://github.com/DRL/blobtools)
- [SnakeMake](http://snakemake.readthedocs.io/en/stable/)
- Additional conda packages not yet listed
