import os

singularity: "docker://genomehubs/blobtoolkit:1.3.3"

include: '../scripts/functions.py'

use_singularity = check_config()

multicore = int(os.getenv('MULTICORE', 1))
maxcore = int(os.getenv('MAXCORE', 1))

asm = config['assembly']['prefix']


rule test:
    """
    Set busco4 output as target of pipeline
    """
    input:
        expand("%s.busco.{lineage}.tsv" % asm, lineage=config['busco']['lineages']),
        expand("%s.busco.{lineage}.txt" % asm, lineage=config['busco']['lineages']),


include: '../rules/run_busco4.smk'
