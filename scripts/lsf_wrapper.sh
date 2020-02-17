#!/bin/bash
#BSUB -n 1
#BSUB -o snakemake.%J.out
#BSUB -e snakemake.%J.err
#BSUB -q long
#BSUB -J insdc-pipeline

# export PIPELINE=/path/to/insdc-pipeline
# export WORKDIR=/path/to/working/directory
# export ASSEMBLY=ABCD01
# export CLUSTER_CONFIG=/path/to/cluster.yaml
# export THREADS=64
# bsub -R "select[mem>10000] rusage[mem=10000]" -M 10000 < $PIPELINE/scripts/lsf_wrapper.sh

# Check the insdc-pipeline directory has been specified
if [ -z "$PIPELINE" ]; then
  echo "ERROR: you must specify PIPELINE=/path/to/insdc-pipeline" >&2
  exit 1
fi

# Check a working directory has been specified
if [ -z "$WORKDIR" ]; then
  echo "ERROR: you must specify WORKDIR=/path/to/working/directory" >&2
  exit 1
fi

# Check an assembly prefix has been specified
if [ -z "$ASSEMBLY" ]; then
  echo "ERROR: you must specify ASSEMBLY=ABCD01" >&2
  exit 1
fi

# Check a cluster.yaml config file has been specified
if [ -z "$CLUSTER_CONFIG" ]; then
  echo "ERROR: you must specify CLUSTER_CONFIG=/path/to/cluster.yaml" >&2
  exit 1
fi

# Check a maximum number of threads has been specified
if [ -z "$THREADS" ]; then
  echo "ERROR: you must specify the maximum number of threads to use, e.g. export THREADS=16" >&2
  exit 1
fi

# Check the specified insdc-pipeline directory contains a Snakefile
if [ ! -s "$PIPELINE/Snakefile" ]; then
  echo "ERROR: $PIPELINE is not a valid insdc-pipeline directory" >&2
  exit 1
fi

# Check the specified working directory exists
if [ ! -d "$WORKDIR" ]; then
  echo "ERROR: $WORKDIR is not a valid working directory" >&2
  exit 1
fi

# Check an assembly config file exists
if [ ! -s "$WORKDIR/$ASSEMBLY.yaml" ]; then
  echo "ERROR: $WORKDIR/$ASSEMBLY.yaml is not a valid configuration file" >&2
  exit 1
fi

# Check a cluster.yaml config file exists
if [ ! -s "$CLUSTER_CONFIG" ]; then
  echo "ERROR: $CLUSTER_CONFIG is not a valid cluster configuration file" >&2
  exit 1
fi

# Load the Conda environment containing snakemake
source $HOME/miniconda3/bin/activate snake_env
EXITCODE=$?
if [ ! $EXITCODE -eq 0 ]; then
  echo "ERROR: Unable to activate conda environment 'snake_env'" >&2
  exit 1
fi

# Check activate binary is available
ACTIVATE=$(which activate)
if [ -z "$ACTIVATE" ]; then
  echo "ERROR: The activate binary could not be found" >&2
  echo "       Try copying the binary into the snake_env environment manually:" >&2
  echo "           cp $HOME/miniconda3/bin/activate $HOME/miniconda3/envs/snake_env/bin/" >&2
  exit 1
fi

# Check the working directory is unlocked in case a previous run failed
if [ -d "$WORKDIR/.snakemake" ]; then
  snakemake --directory $WORKDIR \
            --configfile $WORKDIR/$ASSEMBLY.yaml \
            --unlock \
            -s $PIPELINE/Snakefile
fi

# run snakemake workflow
snakemake -p --use-conda \
             --conda-prefix $HOME/.conda \
             --directory $WORKDIR \
             --configfile $WORKDIR/$ASSEMBLY.yaml \
             --cluster-config $CLUSTER_CONFIG \
             --drmaa " -o {log}.o \
                       -e {log}.e \
                       -R \"select[mem>{cluster.mem}] rusage[mem={cluster.mem}] span[hosts=1]\" \
                       -M {cluster.mem} \
                       -n {cluster.threads} \
                       -q {cluster.queue} " \
             --stats $ASSEMBLY.snakemake.stats \
             --resources btk=1 \
             --latency-wait 300 \
             --rerun-incomplete \
             -j $THREADS \
             -s $PIPELINE/Snakefile && touch $WORKDIR/$ASSEMBLY.complete

