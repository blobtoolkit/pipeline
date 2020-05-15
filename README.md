# insdc-pipeline

_A [Snakemake](http://snakemake.readthedocs.io/en/stable/) pipeline to run analyses on public genome assemblies for visualisation in the [BlobToolKit Viewer](https://blobtoolkit.genomehubs.org/view/)_

[![DOI](https://zenodo.org/badge/125039192.svg)](https://zenodo.org/badge/latestdoi/125039192)

While designed for use on public genome assemblies and read files, this pipeline is also suitable for use with local assembly and read files. If the assembly and read files exist locally in the working directory when the pipeline is run, the remote fetches will be skipped.

More information and tutorials are available at [blobtoolkit.genomehubs.org/pipeline/](https://blobtoolkit.genomehubs.org/pipeline/)

## Overview

The final output of this pipeline is a `BlobDir` dataset directory containing a set of JSON files ready for interactive visualisation with the [BlobToolKit Viewer](https://blobtoolkit.genomehubs.org/view/) (Figure 1).

![Figure 1](https://blobtoolkit.genomehubs.org/wp-content/uploads/2018/08/figure1.png)

The pipeline (Figure 2) is implemented using Snakemake and will fetch all required database and assembly files then run BLAST/Diamond similarity searches and bwa/minimap2 read mapping and for processing with BlobTools.

![Figure 2](https://blobtoolkit.genomehubs.org/wp-content/uploads/2019/11/Figure_2-1.png)

## Installation

Installation instructions for all BlobToolKit components have moved to [blobtoolkit.genomehubs.org/install](https://blobtoolkit.genomehubs.org/install/)

## Usage

Information on pipeline configuration and running has been moved to [blobtoolkit.genomehubs.org/pipeline/pipeline-tutorials](https://blobtoolkit.genomehubs.org/pipeline/pipeline-tutorials/)

## Format specification

The `BlobDir` output of this pipeline contains a set of JSON files and can be validated using the associated JSON-schema. A complete specification and validator for the `BlobDir` format are available at [github.com/blobtoolkit/specification](https://github.com/blobtoolkit/specification).
