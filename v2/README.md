# BTK Pipeline v2

Splits original pipeline into mini-pipelines

## Usage

```
git clone -b develop https://github.com/insdc-pipeline
```

### Generate config files

For public assemblies, the script `generate config.py` will fetch copies of assembly and read files and generate a YAML config file.

```
conda create -y -n btq_env -f envs/btq.yml
conda activate btq_env

DATA_DIR=/path/to/data
DATABASE_DIR=/path/to/databases

./scripts/generate_config.py $ACCESSION --out $DATA_DIR --db $DATABASE_DIR
```

### Run pipelines

```
conda create -y -n btk_env -f envs/btk.yml
conda activate btk_env

TOOL=minimap
THREADS=32

snakemake -p --directory $DATA_DIR/$ACCESSION/$TOOL --configfile $DATA_DIR/GCA_903992545.1/config.yaml -s $TOOL.smk -j $THREADS
```
