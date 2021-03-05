# BTK Pipeline v2

Splits original pipeline into mini-pipelines

## Usage

```
git clone -b develop https://github.com/blobtoolkit/insdc-pipeline
cd insdc-pipeline/v2
```

### Generate config files

For public assemblies, the script `generate config.py` will fetch copies of assembly and read files and generate a YAML config file.

```
conda create -n btq_env -c tolkit tolkein
conda activate btq_env
conda install --clobber -c bioconda entrez-direct
conda install defusedxml

DATA_DIR=/path/to/data
DATABASE_DIR=/path/to/databases

./scripts/generate_config.py $ACCESSION --out $DATA_DIR --db $DATABASE_DIR --coverage 30
```

### Run pipelines

```
conda create -n bts_env -c bioconda snakemake
conda activate bts_env
conda install --clobber -c conda-forge -c bioconda busco=5
conda install --clobber -c bioconda samtools=1.10 minimap2 seqtk
conda activate bts_env

export TOOL=minimap
export THREADS=32
export SNAKE_DIR=`pwd`
export DATA_DIR=/path/to/data
export ACCESSION=GCA_0000000.0
export TOOL=busco

bsub < ./scripts/lsf_wrapper.sh
```
