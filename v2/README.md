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

./scripts/generate_config.py $ACCESSION --out $DATA_DIR --db $DATABASE_DIR
```

### Run pipelines

```
conda env create -n btk_env -f envs/btk.yml
conda activate btk_env

TOOL=minimap
THREADS=32

snakemake -p --directory $DATA_DIR/$ACCESSION/$TOOL --configfile $DATA_DIR/GCA_903992545.1/config.yaml -s $TOOL.smk -j $THREADS
```
