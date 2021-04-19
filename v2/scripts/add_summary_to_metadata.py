#!/usr/bin/env python3

import logging
import os
import re
import subprocess
import traceback
from collections import OrderedDict

import git
import yaml

from functions import read_similarity_settings, reads_by_prefix

logger_config = {
    "level": logging.INFO,
    "format": "%(asctime)s [%(levelname)s] line %(lineno)d %(message)s",
    "filemode": "w",
}
try:
    logger_config.update({"filename": snakemake.log[0]})
except NameError as err:
    pass
logging.basicConfig(**logger_config)
logger = logging.getLogger()

logger.info("Starting script: " + __file__)

try:
    CONFIG = snakemake.config
    ASSEMBLY = snakemake.wildcards.assembly
    OUTFILE = str(snakemake.output)
    GITDIR = (
        git.Repo(
            os.path.dirname(os.path.abspath(__file__)), search_parent_directories=True
        ).working_tree_dir
        + "/.git"
    )
except Exception as err:
    logger.error(traceback.format_exc())
    exit(7)

try:
    meta = {}
    config = CONFIG
    meta["assembly"] = config["assembly"]
    meta["taxon"] = config["taxon"]
    meta["settings"] = config["settings"]
    meta["similarity"] = {
        "diamond_blastx": read_similarity_settings(config, "diamond_blastx"),
        "diamond_blastp": read_similarity_settings(config, "diamond_blastp"),
    }
    meta["reads"] = reads_by_prefix(config)
    for key, value in meta["reads"].items():
        if "url" in value and not isinstance(value["url"], list):
            value["url"] = value["url"].split(";")
        elif "url" not in value:
            value["url"] = []
    p = subprocess.Popen(
        ["git", "--git-dir", GITDIR, "rev-parse", "--short", "HEAD"],
        stdout=subprocess.PIPE,
        encoding="utf-8",
    )
    (output, err) = p.communicate()
    p_status = p.wait()
    meta["settings"]["commit"] = output.rstrip("\n")
    p = subprocess.Popen(
        ["git", "--git-dir", GITDIR, "remote", "-v"],
        stdout=subprocess.PIPE,
        encoding="utf-8",
    )
    (output, err) = p.communicate()
    p_status = p.wait()
    try:
        meta["settings"]["pipeline"] = output.split()[1]
    except IndexError:
        meta["settings"]["pipeline"] = "UNKNOWN"

    yaml.add_representer(
        OrderedDict,
        lambda dumper, data: dumper.represent_mapping(
            "tag:yaml.org,2002:map", data.items()
        ),
    )
    yaml.add_representer(
        tuple,
        lambda dumper, data: dumper.represent_sequence("tag:yaml.org,2002:seq", data),
    )
    yaml.Dumper.ignore_aliases = lambda *args: True
    with open(OUTFILE, "w") as fh:
        fh.write(yaml.dump(meta, indent=1))

except Exception as err:
    logger.error(traceback.format_exc())
    exit(13)
